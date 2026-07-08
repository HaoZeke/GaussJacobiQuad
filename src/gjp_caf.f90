!> Coarray Fortran helpers and multi-rule batching.
!>
!> Single-rule recurrence already partitions Newton roots in gjp_rec (rec_caf).
!> Golub–Welsch and algo665 are dominated by a serial tridiagonal eigensolve, so
!> multi-image speedup for those methods comes from independent *rules*
!> (different α/β or independent problems) partitioned across images.
module gjp_caf
use gjp_types, only: dp
use gjp_rec, only: gauss_jacobi_rec, gauss_jacobi_rec_caf, caf_root_owner, caf_owns_index
use gjp_gw, only: gauss_jacobi_gw
use gjp_algo665, only: gauss_jacobi_algo665
implicit none

private
public :: caf_root_owner, caf_owns_index
public :: gauss_jacobi_batch_caf
public :: gauss_jacobi_one_serial

contains

!> Dispatch one serial (or node-CAF for rec_caf) rule.
subroutine gauss_jacobi_one_serial(npts, alpha, beta, x, wts, method)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    character(len=*), intent(in) :: method

    select case (trim(method))
    case ("rec")
        call gauss_jacobi_rec(npts, alpha, beta, x, wts)
    case ("rec_caf")
        ! Nested CAF node partition if multi-image; still valid inside a batch
        call gauss_jacobi_rec_caf(npts, alpha, beta, x, wts)
    case ("gw", "gw_caf")
        call gauss_jacobi_gw(npts, alpha, beta, x, wts)
    case ("algo665", "algo665_caf")
        call gauss_jacobi_algo665(npts, alpha, beta, x, wts)
    case default
        error stop "gauss_jacobi_one_serial: unknown method"
    end select
end subroutine gauss_jacobi_one_serial

!> Compute nbatch independent rules (same npts, varying α/β), partitioning
!> batch indices across CAF images and gathering with coarrays.
!>
!> Methods: "rec", "gw", "algo665" (serial kernel per owned rule).
!> For "rec_caf", each owned rule uses node-level CAF (usually prefer "rec"
!> inside a large batch so images are not nested-oversubscribed).
subroutine gauss_jacobi_batch_caf(npts, nbatch, alphas, betas, xs, wts, method)
    integer, intent(in) :: npts, nbatch
    real(dp), intent(in) :: alphas(nbatch), betas(nbatch)
    real(dp), intent(out) :: xs(npts, nbatch), wts(npts, nbatch)
    character(len=*), intent(in) :: method

    real(dp), allocatable :: x_caf(:, :)[:], w_caf(:, :)[:]
    real(dp) :: x_local(npts), w_local(npts)
    integer :: me, nimg, i, img
    character(len=:), allocatable :: kernel

    if (nbatch <= 0) error stop "gauss_jacobi_batch_caf: nbatch must be positive"
    if (npts <= 0) error stop "gauss_jacobi_batch_caf: npts must be positive"

    me = this_image()
    nimg = num_images()

    ! Inside a batch, use serial rec so ownership is only over rules
    select case (trim(method))
    case ("rec", "rec_caf")
        kernel = "rec"
    case ("gw", "gw_caf")
        kernel = "gw"
    case ("algo665", "algo665_caf")
        kernel = "algo665"
    case default
        error stop "gauss_jacobi_batch_caf: unknown method"
    end select

    allocate (x_caf(npts, nbatch)[*], w_caf(npts, nbatch)[*])
    x_caf = 0.0_dp
    w_caf = 0.0_dp

    do i = 1, nbatch
        if (.not. caf_owns_index(i, me, nimg)) cycle
        call gauss_jacobi_one_serial(npts, alphas(i), betas(i), x_local, w_local, kernel)
        x_caf(:, i) = x_local
        w_caf(:, i) = w_local
    end do

    sync all

    ! Image 1 gathers each rule from its owner into the local result arrays,
    ! writes the full batch back into its coarray segment, then every image
    ! pulls the complete batch from image 1 (avoids co_broadcast portability issues).
    if (me == 1) then
        do i = 1, nbatch
            img = caf_root_owner(i, nimg)
            if (img == 1) then
                xs(:, i) = x_caf(:, i)
                wts(:, i) = w_caf(:, i)
            else
                xs(:, i) = x_caf(:, i)[img]
                wts(:, i) = w_caf(:, i)[img]
            end if
        end do
        x_caf = xs
        w_caf = wts
    end if
    sync all
    xs = x_caf(:, :)[1]
    wts = w_caf(:, :)[1]

    sync all
    deallocate (x_caf, w_caf)
end subroutine gauss_jacobi_batch_caf

end module gjp_caf

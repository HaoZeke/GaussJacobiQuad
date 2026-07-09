!> Coarray Fortran helpers and multi-rule batching.
!>
!> Single-rule node/index CAF: rec_caf, sturm_caf, glr_caf, bogaert_caf, algo665_dc_caf.
!> Batch CAF partitions independent *rules* (α/β) across images for many kernels.
module gjp_caf
use gjp_types, only: dp
use gjp_rec, only: gauss_jacobi_rec, gauss_jacobi_rec_caf, caf_root_owner, caf_owns_index
use gjp_gw, only: gauss_jacobi_gw
use gjp_algo665, only: gauss_jacobi_algo665, gauss_jacobi_algo665_dc
use gjp_sturm, only: gauss_jacobi_sturm, gauss_jacobi_sturm_caf
use gjp_glr, only: gauss_jacobi_glr, gauss_jacobi_glr_caf
use gjp_bogaert, only: gauss_jacobi_bogaert, gauss_jacobi_bogaert_caf
implicit none

private
public :: caf_root_owner, caf_owns_index
public :: gauss_jacobi_batch_caf
public :: gauss_jacobi_one_serial

contains

!> Dispatch one serial (or single-rule CAF) kernel by name.
subroutine gauss_jacobi_one_serial(npts, alpha, beta, x, wts, method)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    character(len=*), intent(in) :: method

    select case (trim(method))
    case ("rec")
        call gauss_jacobi_rec(npts, alpha, beta, x, wts)
    case ("rec_caf")
        call gauss_jacobi_rec_caf(npts, alpha, beta, x, wts)
    case ("gw", "gw_caf")
        call gauss_jacobi_gw(npts, alpha, beta, x, wts)
    case ("algo665", "algo665_caf")
        call gauss_jacobi_algo665(npts, alpha, beta, x, wts)
    case ("algo665_dc")
        call gauss_jacobi_algo665_dc(npts, alpha, beta, x, wts, use_caf=.false.)
    case ("algo665_dc_caf")
        call gauss_jacobi_algo665_dc(npts, alpha, beta, x, wts, use_caf=.true.)
    case ("sturm")
        call gauss_jacobi_sturm(npts, alpha, beta, x, wts)
    case ("sturm_caf")
        call gauss_jacobi_sturm_caf(npts, alpha, beta, x, wts)
    case ("glr")
        call gauss_jacobi_glr(npts, alpha, beta, x, wts)
    case ("glr_caf")
        call gauss_jacobi_glr_caf(npts, alpha, beta, x, wts)
    case ("bogaert")
        call gauss_jacobi_bogaert(npts, alpha, beta, x, wts)
    case ("bogaert_caf")
        call gauss_jacobi_bogaert_caf(npts, alpha, beta, x, wts)
    case default
        error stop "gauss_jacobi_one_serial: unknown method"
    end select
end subroutine gauss_jacobi_one_serial

!> Compute nbatch independent rules, partitioning batch indices across images.
!> Inside a batch, serial kernels are preferred (avoid nested CAF oversubscribe).
subroutine gauss_jacobi_batch_caf(npts, nbatch, alphas, betas, xs, wts, method)
    integer, intent(in) :: npts, nbatch
    real(dp), intent(in) :: alphas(nbatch), betas(nbatch)
    real(dp), intent(out) :: xs(npts, nbatch), wts(npts, nbatch)
    character(len=*), intent(in) :: method

    real(dp), allocatable :: x_caf(:, :) [:], w_caf(:, :) [:]
    real(dp) :: x_local(npts), w_local(npts)
    integer :: me, nimg, i, img
    character(len=:), allocatable :: kernel

    if (nbatch <= 0) error stop "gauss_jacobi_batch_caf: nbatch must be positive"
    if (npts <= 0) error stop "gauss_jacobi_batch_caf: npts must be positive"

    me = this_image()
    nimg = num_images()

    select case (trim(method))
    case ("rec", "rec_caf")
        kernel = "rec"
    case ("gw", "gw_caf")
        kernel = "gw"
    case ("algo665", "algo665_caf")
        kernel = "algo665"
    case ("algo665_dc", "algo665_dc_caf")
        kernel = "algo665_dc"
    case ("sturm", "sturm_caf")
        kernel = "sturm"
    case ("glr", "glr_caf")
        kernel = "glr"
    case ("bogaert", "bogaert_caf")
        kernel = "bogaert"
    case default
        error stop "gauss_jacobi_batch_caf: unknown method"
    end select

    allocate (x_caf(npts, nbatch) [*], w_caf(npts, nbatch) [*])
    x_caf = 0.0_dp
    w_caf = 0.0_dp

    do i = 1, nbatch
        if (.not. caf_owns_index(i, me, nimg)) cycle
        call gauss_jacobi_one_serial(npts, alphas(i), betas(i), x_local, w_local, kernel)
        x_caf(:, i) = x_local
        w_caf(:, i) = w_local
    end do

    sync all

    if (me == 1) then
        do i = 1, nbatch
            img = caf_root_owner(i, nimg)
            if (img == 1) then
                xs(:, i) = x_caf(:, i)
                wts(:, i) = w_caf(:, i)
            else
                xs(:, i) = x_caf(:, i) [img]
                wts(:, i) = w_caf(:, i) [img]
            end if
        end do
        x_caf = xs
        w_caf = wts
    end if
    sync all
    xs = x_caf(:, :) [1]
    wts = w_caf(:, :) [1]

    sync all
    deallocate (x_caf, w_caf)
end subroutine gauss_jacobi_batch_caf

end module gjp_caf

!> Batch CAF vs serial for rec, gw, and algo665 (shipped batch path).
program test_caf_batch
use GaussJacobiQuad, only: gauss_jacobi, gauss_jacobi_batch_caf
use gjp_types, only: dp
implicit none

integer, parameter :: npts = 24, nbatch = 8
real(dp) :: alphas(nbatch), betas(nbatch)
real(dp) :: xs(npts, nbatch), ws(npts, nbatch)
real(dp) :: xs_ref(npts, nbatch), ws_ref(npts, nbatch)
real(dp) :: max_dx, max_dw, tol
integer :: ib, me, nimg, imethod
character(len=16) :: methods(3)
character(len=:), allocatable :: method
logical :: ok

me = this_image()
nimg = num_images()
ok = .true.
tol = 1.0e-12_dp
methods = [character(len=16) :: "rec", "gw", "algo665"]

do ib = 1, nbatch
    alphas(ib) = 0.5_dp
    betas(ib) = 0.0_dp + 0.05_dp * real(ib - 1, dp)
end do

do imethod = 1, 3
    method = trim(methods(imethod)) ! deferred length
    do ib = 1, nbatch
        call gauss_jacobi(npts, alphas(ib), betas(ib), xs_ref(:, ib), ws_ref(:, ib), method)
    end do
    call gauss_jacobi_batch_caf(npts, nbatch, alphas, betas, xs, ws, method)
    max_dx = maxval(abs(xs - xs_ref))
    max_dw = maxval(abs(ws - ws_ref))
    if (me == 1) then
        print '(A,A,A,I0,A,ES12.5,A,ES12.5)', &
            "batch method=", method, " nimg=", nimg, &
            " max|dx|=", max_dx, " max|dw|=", max_dw
    end if
    if (max_dx > tol .or. max_dw > tol) ok = .false.
end do

sync all
if (.not. ok) error stop "test_caf_batch failed"
if (me == 1) print*,"test_caf_batch PASSED (rec/gw/algo665 batch CAF)"
end program test_caf_batch

!> Compare Cuppen+imtqlx DC (and CAF secular) against classic algo665/imtqlx.
program test_imtqlx_dc
use GaussJacobiQuad, only: gauss_jacobi
use gjp_types, only: dp
implicit none

integer, parameter :: n = 40
real(dp), parameter :: alpha = 0.5_dp, beta = 1.5_dp
real(dp) :: x_ref(n), w_ref(n), x_dc(n), w_dc(n), x_caf(n), w_caf(n)
real(dp) :: max_dx, max_dw, tol
integer :: me
character(len=:), allocatable :: m_ref, m_dc, m_caf

me = this_image()
tol = 1.0e-12_dp
m_ref = "algo665"
m_dc = "algo665_dc"
m_caf = "algo665_dc_caf"

call gauss_jacobi(n, alpha, beta, x_ref, w_ref, m_ref)
call gauss_jacobi(n, alpha, beta, x_dc, w_dc, m_dc)
call gauss_jacobi(n, alpha, beta, x_caf, w_caf, m_caf)

! Sort-independent compare: both should be sorted ascending
max_dx = maxval(abs(x_dc - x_ref))
max_dw = maxval(abs(w_dc - w_ref))
if (me == 1) then
    print '(A,ES12.5,A,ES12.5)', "algo665 vs algo665_dc: max|dx|=", max_dx, " max|dw|=", max_dw
end if
if (max_dx > tol .or. max_dw > tol) then
    if (me == 1) print *, "FAIL dc vs algo665"
    error stop "imtqlx_dc disagree"
end if

max_dx = maxval(abs(x_caf - x_ref))
max_dw = maxval(abs(w_caf - w_ref))
if (me == 1) then
    print '(A,ES12.5,A,ES12.5)', "algo665 vs algo665_dc_caf: max|dx|=", max_dx, " max|dw|=", max_dw
    print '(A,I0)', "test_imtqlx_dc PASSED nimg=", num_images()
end if
if (max_dx > tol .or. max_dw > tol) error stop "imtqlx_dc_caf disagree"
end program test_imtqlx_dc

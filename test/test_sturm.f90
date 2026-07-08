!> Sturm (+ CAF) vs algo665 on the same Golub-Welsch Jacobi matrix.
program test_sturm
use GaussJacobiQuad, only: gauss_jacobi
use gjp_types, only: dp
implicit none

integer, parameter :: n = 48
real(dp), parameter :: alpha = 0.5_dp, beta = 2.0_dp
real(dp) :: x_ref(n), w_ref(n), x_s(n), w_s(n), x_c(n), w_c(n)
real(dp) :: max_dx, max_dw, tol
integer :: me
character(len=:), allocatable :: m_ref, m_s, m_c

me = this_image()
tol = 1.0e-8_dp
m_ref = "algo665"
m_s = "sturm"
m_c = "sturm_caf"

call gauss_jacobi(n, alpha, beta, x_ref, w_ref, m_ref)
call gauss_jacobi(n, alpha, beta, x_s, w_s, m_s)
call gauss_jacobi(n, alpha, beta, x_c, w_c, m_c)

max_dx = maxval(abs(x_s - x_ref))
max_dw = maxval(abs(w_s - w_ref))
if (me == 1) print '(A,ES12.5,A,ES12.5)', "sturm vs algo665 max|dx|=", max_dx, " max|dw|=", max_dw
if (max_dx > tol .or. max_dw > tol) error stop "sturm disagree"

max_dx = maxval(abs(x_c - x_ref))
max_dw = maxval(abs(w_c - w_ref))
if (me == 1) then
    print '(A,ES12.5,A,ES12.5)', "sturm_caf vs algo665 max|dx|=", max_dx, " max|dw|=", max_dw
    print '(A,I0)', "test_sturm PASSED nimg=", num_images()
end if
if (max_dx > tol .or. max_dw > tol) error stop "sturm_caf disagree"
end program test_sturm

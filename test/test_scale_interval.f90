!> Integrate a polynomial on a non-symmetric interval with the scaled rule.
program test_scale_interval
use GaussJacobiQuad, only: gauss_jacobi_rule, gauss_jacobi_scale
use gjp_types, only: dp
implicit none

integer, parameter :: n = 8
real(dp), parameter :: alpha = 1.0_dp, beta = 2.0_dp
real(dp), parameter :: a = 0.0_dp, b = 3.0_dp
real(dp) :: x(n), w(n), quad, exact, err, wsum
integer :: i

call gauss_jacobi_rule(n, alpha, beta, x, w, "rec")
call gauss_jacobi_scale(n, alpha, beta, a, b, x, w)

quad = 0.0_dp
do i = 1, n
    quad = quad + w(i) * x(i)**3
end do
! ∫_0^3 (3-t)^1 t^2 t^3 dt = ∫_0^3 (3-t) t^5 dt = 729/14
exact = 729.0_dp / 14.0_dp
err = abs(quad - exact)
! f = 1 moment: (b-a)^(alpha+beta+1) * B(beta+1, alpha+1) = 81 * B(3,2) = 27/4
wsum = sum(w)
print '(A,ES16.8,A,ES16.8)', "scale [0,3] poly err=", err, " moment err=", abs(wsum - 27.0_dp / 4.0_dp)
if (err > 1.0e-12_dp) error stop "scaled polynomial integral failed"
if (abs(wsum - 27.0_dp / 4.0_dp) > 1.0e-12_dp) error stop "scaled weight sum failed"
if (any(x <= a) .or. any(x >= b)) error stop "scaled nodes outside (a,b)"
if (any(w <= 0.0_dp)) error stop "scaled weights not positive"
end program test_scale_interval

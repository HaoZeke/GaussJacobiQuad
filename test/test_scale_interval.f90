!> Integrate a polynomial on a non-symmetric interval with the scaled rule.
program test_scale_interval
use GaussJacobiQuad, only: gauss_jacobi_rule, gauss_jacobi_scale
use gjp_types, only: dp
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
implicit none

integer, parameter :: n = 8
integer, parameter :: n4 = 4
real(dp), parameter :: alpha = 1.0_dp, beta = 2.0_dp
real(dp), parameter :: a = 0.0_dp, b = 3.0_dp
real(dp) :: x(n), w(n), quad, exact, err, wsum
real(dp) :: x4(n4), w4(n4), y4(n4), z4(n4), lo, hi
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

! Huge finite endpoints: finite nodes and weights inside (a, b).
! Status is not returned here; a rejected scale error-stops.
call gauss_jacobi_rule(n4, 0.0_dp, 0.0_dp, x4, w4, "rec")
y4 = x4
z4 = w4
lo = 1.0e308_dp
hi = 1.7e308_dp
call gauss_jacobi_scale(n4, 0.0_dp, 0.0_dp, lo, hi, x4, w4)
if (.not. all(ieee_is_finite(x4)) .or. .not. all(ieee_is_finite(w4))) &
    error stop "non-finite values on [1e308, 1.7e308]"
if (any(x4 <= lo) .or. any(x4 >= hi)) error stop "nodes outside [1e308, 1.7e308]"
lo = -1.0e308_dp
hi = 1.0e308_dp
call gauss_jacobi_scale(n4, 0.0_dp, 0.0_dp, lo, hi, y4, z4)
if (.not. all(ieee_is_finite(y4)) .or. .not. all(ieee_is_finite(z4))) &
    error stop "non-finite values on [-1e308, 1e308]"
if (any(y4 <= lo) .or. any(y4 >= hi)) error stop "nodes outside [-1e308, 1e308]"
end program test_scale_interval

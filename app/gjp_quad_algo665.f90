program gjp_quad_algo665

use GaussJacobiQuad, only: gauss_jacobi
use gjp_types, only: dp
implicit none

integer :: n_points
real(dp) :: alpha, beta
real(dp), dimension(:), allocatable :: x, w
character(len=128) :: arg
character(len=:), allocatable :: method
integer :: idx

if (command_argument_count() /= 3) then
    print*,"Usage: ./gjp_quad_algo665 <n_points> <alpha> <beta>"
    print*,"  n_points: Number of quadrature points (integer)"
    print*,"  alpha: Parameter alpha for Gauss-Jacobi quadrature (must be > -1)"
    print*,"  beta: Parameter beta for Gauss-Jacobi quadrature (must be > -1)"
    error stop "Must supply 3 arguments"
end if

call get_command_argument(1, arg)
read (arg, '(i4)') n_points
allocate (x(n_points), w(n_points))

call get_command_argument(2, arg)
if (index(arg, '.') == 0) then
    print*,"alpha must include a decimal point"
    error stop
end if
read (arg, *) alpha

call get_command_argument(3, arg)
if (index(arg, '.') == 0) then
    print*,"beta must include a decimal point"
    error stop
end if
read (arg, *) beta

method = "algo665"
call gauss_jacobi(n_points, alpha, beta, x, w, method)

do idx = 1, n_points
    print '(1X, A, 1P, E24.17, 2X, A, 1P, E23.17)', 'Root: ', x(idx), 'Weight: ', w(idx)
end do

deallocate (x, w)

end program gjp_quad_algo665

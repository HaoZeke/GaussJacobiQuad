program gjp_quad_rec

use GaussJacobiQuad, only: gauss_jacobi_rec
use gjp_types, only: dp
implicit none

integer :: n_points
real(dp) :: alpha, beta
real(dp), dimension(:), allocatable :: x, w
character(len=128) :: arg
integer :: i

if (command_argument_count() /= 3) then
    print*,"./gjp_quad <n_points> <alpha> <beta>"
    error stop "Must supply 3 arguments"
end if

call get_command_argument(1, arg)
read (arg, '(i4)') n_points
allocate (x(n_points), w(n_points))

call get_command_argument(2, arg)
read (arg, '(f8.3)') alpha

call get_command_argument(3, arg)
read (arg, '(f8.3)') beta

call gauss_jacobi_rec(n_points, alpha, beta, x, w)

do i = 1, n_points
    print '(1X, A, 1P, E24.17, 2X, A, 1P, E23.17)', 'Root: ', x(i), 'Weight: ', w(i)
end do

deallocate (x, w)

end program gjp_quad_rec

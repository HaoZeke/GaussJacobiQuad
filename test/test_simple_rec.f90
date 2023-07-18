program test_simple_rec
use GaussJacobiQuad
implicit none
integer, parameter :: n = 5
real(dp), parameter :: a = 0.0_dp, b = 12.0_dp
real(dp), dimension(n) :: x, w
integer :: i
open (unit=10, file='gauss_jacobi_output.txt', status='unknown')

call gauss_jacobi(n, a, b, x, w)

do i = 1, n
    write (10, *) 'Root: ', x(i), ' Weight: ', w(i)
    print '(1X, A, 1P, E24.17, 2X, A, 1P, E23.17)', 'Root: ', x(i), 'Weight: ', w(i)
end do

close (10)
end program test_simple_rec

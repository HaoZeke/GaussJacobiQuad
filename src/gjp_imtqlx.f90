module gjp_imtqlx
use gjp_lapack, only: drot, drotg
use gjp_types, only: dp
implicit none

contains

!> @brief Implicitly-shifted Modified QL factorization algorithm for symmetric tridiagonal matrices.
!>
!> This subroutine diagonalizes a real, symmetric, tridiagonal matrix using the implicitly-shifted Modified QL factorization algorithm.
!> The eigenvalues are overwritten on the diagonal of the matrix, while the eigenvectors are stored in the array `z`. This is an internal
!> function typically used in the context of computing Gauss-Jacobi quadrature points and weights.
!>
!> This function is an adaptation of the Implicitly-shifted QL algorithm. It is specifically optimized for tridiagonal matrices,
!> which is the form that Jacobi matrices often take. This is a modernized version of the algorithm presented in Algorithm 655
!> which in turn is a slightly modified variant of an EISPACK routine.
!>
!> @param[in] n The size of the matrix, i.e., number of rows or columns.
!> @param[inout] d On entry, the diagonal elements of the tridiagonal matrix. On exit, the eigenvalues in ascending order.
!> @param[inout] e On entry, the off-diagonal elements of the tridiagonal matrix. On exit, the values are overwritten.
!> @param[inout] z On entry, the initial guesses for the eigenvectors. On exit, the computed eigenvectors.
!>
!> @note This subroutine has an iteration limit (`itn`) set to 30. An error is reported if the limit is exceeded.
!>
!> @warning Make sure to allocate sufficient space for the arrays `d`, `e`, and `z` before calling this function.
subroutine imtqlx(n, d, e, z)
    implicit none

    ! Inputs
    integer, intent(in) :: n
    real(dp), intent(inout) :: d(n)
    real(dp), intent(inout) :: e(n)
    real(dp), intent(inout) :: z(n)

    ! Local variables
    real(dp) :: prec, p, g, r, s, f, b, c
    integer :: l, m, ii, i, j, k
    integer, parameter :: max_iter = 30
    real(dp), allocatable :: Bjk(:, :)
    integer :: nu, JK

    ! Initialize precision
    prec = epsilon(prec)

    ! Allocate and initialize Bjk
    allocate (Bjk(n, 2 * max_iter))
    Bjk = 0.0_dp

    ! Set the last off-diagonal element to zero
    e(n) = 0.0_dp

    ! Main loop
    do l = 1, n
        j = 0
        do while (j < max_iter)
            do m = l, n
                if (m == n) exit
                if (abs(e(m)) <= prec * (abs(d(m)) + abs(d(m + 1)))) exit
            end do

            p = d(l)
            if (m == l) exit

            if (j > max_iter) then
                print*," "
                print*,"IMTQLX - Fatal error!"
                print*,"Iteration limit exceeded."
                print*,"J = ", j
                print*,"L = ", l
                print*,"M = ", m
                print*,"N = ", n
                error stop "Terminating due to iteration limit exceeded."
            end if

            j = j + 1
            g = (d(l + 1) - p) / (2.0_dp * e(l))
            r = sqrt(g * g + 1.0_dp)
            g = d(m) - p + e(l) / (g + sign(r, g))
            s = 1.0_dp
            c = 1.0_dp
            p = 0.0_dp
            do ii = 1, m - l
                i = m - ii
                f = s * e(i)
                b = c * e(i)

                if (abs(g) <= abs(f)) then
                    c = g / f
                    r = sqrt(c * c + 1.0_dp)
                    e(i + 1) = f * r
                    s = 1.0_dp / r
                    c = c * s
                else
                    s = f / g
                    r = sqrt(s * s + 1.0_dp)
                    e(i + 1) = g * r
                    c = 1.0_dp / r
                    s = s * c
                end if

                g = d(i + 1) - p
                r = (d(i) - g) * s + 2.0_dp * c * b
                p = s * r
                d(i + 1) = g + p
                g = c * r - b
                f = z(i + 1)
                z(i + 1) = s * z(i) + c * f
                z(i) = c * z(i) - s * f
            end do

            d(l) = d(l) - p
            e(l) = g
            e(m) = 0.0_dp
        end do
    end do

    ! Clean up
    deallocate (Bjk)

    ! Sorting
    ! do ii = 2, n
    !     i = ii - 1
    !     k = i
    !     p = d(i)
    !     do j = ii, n
    !         if (d(j) < p) then
    !             k = j
    !             p = d(j)
    !         end if
    !     end do
    !     if (k /= i) then
    !         d(k) = d(i)
    !         d(i) = p
    !         p = z(i)
    !         z(i) = z(k)
    !         z(k) = p
    !     end if
    ! end do
    call dsort2a(n, d, z)
end subroutine

! Sorts x and performs the same swaps on y
subroutine dsort2a(n, x, y)
    integer, intent(in) :: n
    real(dp), intent(inout) :: x(n), y(n)
    integer :: i, j, min_idx
    real(dp) :: temp

    do i = 1, n - 1
        min_idx = i
        do j = i + 1, n
            if (x(j) < x(min_idx)) min_idx = j
        end do
        if (min_idx /= i) then
            ! Swap x(i) and x(min_idx)
            temp = x(i)
            x(i) = x(min_idx)
            x(min_idx) = temp

            ! Swap y(i) and y(min_idx)
            temp = y(i)
            y(i) = y(min_idx)
            y(min_idx) = temp
        end if
    end do
end subroutine dsort2a

end module gjp_imtqlx

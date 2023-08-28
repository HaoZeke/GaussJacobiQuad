! BEGIN_HEADER
! -----------------------------------------------------------------------------
! Gauss-Jacobi Quadrature Implementation
! Authors: Rohit Goswami <rgoswami[at]ieee.org>
! Source: GaussJacobiQuad Library
! License: MIT
! GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
! Date: 2023-08-28
! Commit: c442f77
! -----------------------------------------------------------------------------
! This code is part of the GaussJacobiQuad library, providing an efficient
! implementation for Gauss-Jacobi quadrature nodes and weights computation.
! -----------------------------------------------------------------------------
! To cite this software:
! Rohit Goswami (2023). HaoZeke/GaussJacobiQuad: v0.1.0.
! Zenodo: https://doi.org/10.5281/ZENODO.8285112
! ---------------------------------------------------------------------
! END_HEADER
module gjp_imtqlx
use gjp_types, only: dp
implicit none

contains

!> @brief Implicitly-shifted Modified QL factorization algorithm for symmetric tridiagonal matrices.
!>
!> This subroutine diagonalizes a real, symmetric, tridiagonal matrix using the implicitly-shifted Modified QL factorization algorithm.
!> The eigenvalues are overwritten on the diagonal of the matrix, while the eigenvectors are stored in the array `sol_vec`. This is an internal
!> function typically used in the context of computing Gauss-Jacobi quadrature points and weights.
!>
!> This function is an adaptation of the Implicitly-shifted QL algorithm. It is specifically optimized for tridiagonal matrices,
!> which is the form that Jacobi matrices often take. This is a modernized version of the algorithm presented in Algorithm 655
!> which in turn is a slightly modified variant of an EISPACK routine.
!>
!> @param[in] mat_size The size of the matrix, i.e., number of rows or columns.
!> @param[inout] diag On entry, the diagonal elements of the tridiagonal matrix. On exit, the eigenvalues in ascending order.
!> @param[inout] off_diag On entry, the off-diagonal elements of the tridiagonal matrix. On exit, the values are overwritten.
!> @param[inout] sol_vec On entry, the initial guesses for the eigenvectors. On exit, the computed eigenvectors.
!>
!> @note This subroutine has an iteration limit (`max_iter`) set to 30. An error is reported if the limit is exceeded.
!>
!> @warning Make sure to allocate sufficient space for the arrays `diag`, `off_diag`, and `sol_vec` before calling this function.
subroutine imtqlx(mat_size, diag, off_diag, sol_vec)
    implicit none

    ! Inputs
    integer, intent(in) :: mat_size
    real(dp), intent(inout) :: diag(mat_size)
    real(dp), intent(inout) :: off_diag(mat_size - 1)
    real(dp), intent(inout) :: sol_vec(mat_size)

    ! Local variables
    real(dp) :: precision, pivot_val, g_val, rot_val, scale_val, f_val, b_val, cos_val
    integer :: lower_bound, upper_bound, inner_i, i, iter_count
    integer, parameter :: max_iter = 30

    ! Initialize machine precision
    precision = epsilon(precision)

    ! Set the last off-diagonal element to zero
    off_diag(mat_size - 1) = 0.0_dp

    ! Main loop over all diagonal elements (Eigenvalues)
    do lower_bound = 1, mat_size
        iter_count = 0 ! Initialize iteration counter for the inner loop

        ! Inner loop for convergence criteria (max iterations)
        do while (iter_count < max_iter)
            ! Loop to find upper bound for the current lower bound
            do upper_bound = lower_bound, mat_size
                if (upper_bound == mat_size) exit
                if (abs(off_diag(upper_bound)) <= precision * (abs(diag(upper_bound)) + abs(diag(upper_bound + 1)))) exit
            end do

            ! Check for convergence
            pivot_val = diag(lower_bound)
            if (upper_bound == lower_bound) exit ! Converged

            ! Check for iteration limit
            if (iter_count > max_iter) then
                print*," "
                print*,"IMTQLX - Fatal error!"
                print*,"Iteration limit exceeded."
                stop "Terminating due to iteration limit exceeded."
            end if

            iter_count = iter_count + 1 ! Update iteration count

            ! Calculation of g_val based on pivot and off-diagonal elements
            g_val = (diag(lower_bound + 1) - pivot_val) / (2.0_dp * off_diag(lower_bound))
            rot_val = sqrt(g_val * g_val + 1.0_dp)
            g_val = diag(upper_bound) - pivot_val + off_diag(lower_bound) / (g_val + sign(rot_val, g_val))

            ! Initialize rotation parameters
            scale_val = 1.0_dp
            cos_val = 1.0_dp
            pivot_val = 0.0_dp

            ! Loop for the implicit QR factorization
            do inner_i = 1, upper_bound - lower_bound
                i = upper_bound - inner_i
                f_val = scale_val * off_diag(i)
                b_val = cos_val * off_diag(i)

                ! Update by Givens rotation
                if (abs(g_val) <= abs(f_val)) then
                    cos_val = g_val / f_val
                    rot_val = sqrt(cos_val * cos_val + 1.0_dp)
                    off_diag(i + 1) = f_val * rot_val
                    scale_val = 1.0_dp / rot_val
                    cos_val = cos_val * scale_val
                else
                    scale_val = f_val / g_val
                    rot_val = sqrt(scale_val * scale_val + 1.0_dp)
                    off_diag(i + 1) = g_val * rot_val
                    cos_val = 1.0_dp / rot_val
                    scale_val = scale_val * cos_val
                end if

                ! Update diagonal and off-diagonal elements
                g_val = diag(i + 1) - pivot_val
                rot_val = (diag(i) - g_val) * scale_val + 2.0_dp * cos_val * b_val
                pivot_val = scale_val * rot_val
                diag(i + 1) = g_val + pivot_val
                g_val = cos_val * rot_val - b_val

                ! Update solution vector
                f_val = sol_vec(i + 1)
                sol_vec(i + 1) = scale_val * sol_vec(i) + cos_val * f_val
                sol_vec(i) = cos_val * sol_vec(i) - scale_val * f_val
            end do

            ! Update diagonal and off-diagonal elements after QR step
            diag(lower_bound) = diag(lower_bound) - pivot_val
            off_diag(lower_bound) = g_val
            off_diag(upper_bound) = 0.0_dp
        end do
    end do

    ! Sort the solution after convergence
    call dsort2a(mat_size, diag, sol_vec)
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

!> @brief This module provides utility functions for computing Jacobi matrices and zeroth moments.
!!
!! @details This module contains essential functions for working with Jacobi polynomials. These
!! functions are particularly useful for various numerical methods that utilize Jacobi polynomials.
!! The module provides a function for generating a Jacobi matrix and another for calculating the
!! zeroth moment of Jacobi polynomials.
!!
!! Functions included are:
!! - jacobi_matrix: Generates a Jacobi matrix based on given alpha and beta parameters.
!! - jacobi_zeroeth_moment: Calculates the zeroth moment of Jacobi polynomials.
module gjp_common
use gjp_types, only: dp, gjp_sparse_matrix
implicit none

contains
!> @brief Computes the Jacobi matrix for given parameters.
!>
!> The Jacobi matrix is computed as:
!> \[
!> J_{i,j} = \begin{cases}
!>   \alpha & \text{if } i = j \\
!>   \beta  & \text{if } |i-j| = 1 \\
!>   0     & \text{otherwise}
!> \end{cases}
!> \]
!>
!> @param[in] n Size of the matrix, number of points
!> @param[in] alpha parameter for Jacobi polynomials
!> @param[in] beta parameter for Jacobi polynomials
!> @return A gjp_sparse_matrix representing the Jacobi matrix
function jacobi_matrix(n, alpha, beta) result(jacmat)
    integer, intent(in) :: n ! Size of the matrix, number of points
    real(dp), intent(in) :: alpha, beta
    type(gjp_sparse_matrix) :: jacmat
    integer :: idx
    real(dp) :: ab, abi, a2b2

    allocate (jacmat%diagonal(n))
    allocate (jacmat%off_diagonal(n - 1))

    ab = alpha + beta
    abi = 2.0_dp + ab
    jacmat%diagonal(1) = (beta - alpha) / abi
    jacmat%off_diagonal(1) = 4.0_dp * (1.0_dp + alpha) * (1.0_dp + beta) &
                             / ((abi + 1.0_dp) * abi * abi)
    a2b2 = beta * beta - alpha * alpha
    do idx = 2, n
        abi = 2.0_dp * idx + ab
        jacmat%diagonal(idx) = a2b2 / ((abi - 2.0_dp) * abi)
        abi = abi**2
        if (idx < n) then
            jacmat%off_diagonal(idx) = 4.0_dp * idx * (idx + alpha) * (idx + beta) &
                                       * (idx + ab) / ((abi - 1.0_dp) * abi)
        end if
    end do
    jacmat%off_diagonal(1:n - 1) = sqrt(jacmat%off_diagonal(1:n - 1))
end function jacobi_matrix

!> @brief Computes the zeroth moment for Jacobi polynomials.
!>
!> The zeroth moment is computed using the formula:
!> \[
!> \text{zmom} = 2^{(\alpha + \beta + 1)} \frac{\Gamma(\alpha + 1) \Gamma(\beta + 1)}{\Gamma(2 + \alpha + \beta)}
!> \]
!> Where \(\Gamma\) is the gamma function.
!>
!> @param[in] alpha parameter for Jacobi polynomials
!> @param[in] beta parameter for Jacobi polynomials
!> @return The zeroth moment value
!> @note The zeroth moment should always be positive
function jacobi_zeroeth_moment(alpha, beta) result(zmom)
    real(dp), intent(in) :: alpha, beta
    real(dp) :: zmom
    real(dp) :: ab, abi

    ab = alpha + beta
    abi = 2.0_dp + ab

    zmom = 2.0_dp**(alpha + beta + 1.0_dp) * exp(log_gamma(alpha + 1.0_dp) &
                                                 + log_gamma(beta + 1.0_dp) - log_gamma(abi))

    if (zmom <= 0.0_dp) then
        error stop "Zeroth moment is not positive but should be"
    end if
end function jacobi_zeroeth_moment

end module gjp_common

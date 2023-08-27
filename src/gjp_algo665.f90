!> @brief This module provides routines for numerical integration using Gauss-Jacobi quadrature.
!! @details The implementation is based on the GW method as implemented in the GaussJacobiQuad
!! with a modernized version of the implicit QL in Algorithm 655 from ACM Collected Algorithms.
!!
!! References:
!!  - Elhay, S., & Kautsky, J. (1987). Algorithm 655: IQPACK: FORTRAN Subroutines for the Weights
!!    of Interpolatory Quadratures. ACM Transactions on Mathematical Software, 13(4), 399-415.
!!    DOI: 10.1145/35078.214351
!!  - Martin, C., & Wilkinson, J. H. (1968). The implicit QL algorithm. Numerische Mathematik, 12(5), 377-383.
!!    DOI: 10.1007/BF02165404module gjp_algo665
module gjp_algo665
use gjp_types, only: dp, gjp_sparse_matrix
use gjp_imtqlx, only: imtqlx
implicit none
contains

!> @brief Computes the zeros and weights for Gauss-Jacobi quadrature.
!>
!> This subroutine computes the zeros (`x`) and weights (`w`) for Gauss-Jacobi quadrature
!> by diagonalizing the Jacobi matrix. The solution involves finding the eigenvalues of the Jacobi matrix,
!> which are the roots of the Jacobi polynomial. Since the Jacobi matrix is a symmetric tridiagonal matrix,
!> the LAPACK DSTEQR routine is used, utilizing its specific form for tridiagonal matrices.
!>
!> The Jacobi matrix is represented as:
!> \[
!> J = \begin{bmatrix}
!>   \alpha & \beta  & 0      & \cdots & 0      \\
!>   \beta  & \alpha & \beta  & \cdots & 0      \\
!>   \vdots & \vdots & \ddots & \ddots & \vdots \\
!>   0      & 0      & \cdots & \beta  & \alpha
!> \end{bmatrix}
!> \]
!>
!> `Z` is initialized as the identity matrix, and the eigenvectors are used to compute the wts.
!>
!> @param[in] npts Number of x
!> @param[in] alpha parameter for Jacobi polynomials
!> @param[in] beta parameter for Jacobi polynomials
!> @param[out] x Zeros of Jacobi polynomials
!> @param[out] wts weights for Gauss-Jacobi quadrature
subroutine gauss_jacobi_algo665(npts, alpha, beta, x, wts)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    real(dp) :: zeroeth_moment
    type(gjp_sparse_matrix) :: jacobi_mat
    real(dp) :: diagonal_elements(npts), &
                off_diagonal_elements(npts - 1), &
                eigenvectors(npts, npts), &
                wodpspace(2 * npts - 2)
    integer :: computation_info, idx

    jacobi_mat = jacobi_matrix(npts, alpha, beta)
    zeroeth_moment = jacobi_zeroeth_moment(alpha, beta)

    ! Extract diagonal and off-diagonal elements
    diagonal_elements = jacobi_mat%diagonal(1:npts)
    off_diagonal_elements = jacobi_mat%off_diagonal(1:npts - 1)

    ! Initialize weights and knot points
    wts = 0.0_dp
    x = diagonal_elements
    wts(1) = sqrt(zeroeth_moment)

    ! Diagonalize the Jacobi matrix using the modified implicitly shifted QL method.
    call imtqlx(npts, x, off_diagonal_elements, wts)

    ! The weights are related to the squares of the first components of the
    ! eigenvectors
    wts = wts**2

end subroutine gauss_jacobi_algo665

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

end module gjp_algo665

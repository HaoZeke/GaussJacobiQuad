! BEGIN_HEADER
! -----------------------------------------------------------------------------
! Gauss-Jacobi Quadrature Implementation
! Authors: Rohit Goswami <rgoswami[at]ieee.org>
! Source: GaussJacobiQuad Library
! License: MIT
! GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
! Date: 2023-08-28
! Commit: 69e8946
! -----------------------------------------------------------------------------
! This code is part of the GaussJacobiQuad library, providing an efficient
! implementation for Gauss-Jacobi quadrature nodes and weights computation.
! -----------------------------------------------------------------------------
! END_HEADER
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
use gjp_common, only: jacobi_matrix, jacobi_zeroeth_moment
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
                off_diagonal_elements(npts - 1)

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

    ! The weights are related to the squares of the the eigenvectors
    wts = wts**2

end subroutine gauss_jacobi_algo665

end module gjp_algo665

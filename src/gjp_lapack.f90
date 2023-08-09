! BEGIN_HEADER
! -----------------------------------------------------------------------------
! Gauss-Jacobi Quadrature Implementation
! Authors: Rohit Goswami <rgoswami[at]ieee.org>
! Source: GaussJacobiQuad Library
! License: MIT
! GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
! Date: 2023-08-09
! Commit: b04e1b3
! -----------------------------------------------------------------------------
! This code is part of the GaussJacobiQuad library, providing an efficient
! implementation for Gauss-Jacobi quadrature nodes and weights computation.
! -----------------------------------------------------------------------------
! END_HEADER
!> @brief LAPACK Interface for Gauss Jacobi Polynomials
!> @details This module provides an interface to LAPACK routines used with Generalized Jacobi Polynomials.
module gjp_lapack
implicit none

integer, parameter :: dp = kind(0.d0)

!> @brief Compute the eigenvalues and eigenvectors of a real symmetric tridiagonal matrix
!> @details The DSTEQR subroutine computes all eigenvalues and, optionally, eigenvectors of a
!> real, symmetric, tridiagonal matrix. It uses the implicit QL or QR method to do so.
!>
!> More information on DSTEQR can be found in the [LAPACK documentation](https://netlib.org/lapack/explore-html/d2/d24/group__aux_o_t_h_e_rcomputational_ga47fe470e7a882c58d4cc49e2c6cb7c70.html).
!>
!> @param COMPZ Specifies computation of eigenvectors ('V') or not ('N')
!> @param N Order of the matrix
!> @param D On entry, the diagonal elements. On exit, the eigenvalues in ascending order
!> @param E Off-diagonal elements on entry; destroyed on exit
!> @param Z On entry, if COMPZ = 'V', then Z is an N-by-N matrix. On exit, if COMPZ = 'V', Z contains the orthogonal matrix of eigenvectors
!> @param LDZ Leading dimension of Z as declared in the calling subroutine
!> @param WORK Work array, dimension (max(1, 2*N-2)) if COMPZ = 'V', (max(1, 1)) otherwise
!> @param INFO INFO = 0: successful exit. INFO = K: the algorithm failed to find all the eigenvalues, no eigenvectors or eigenvectors of order K have been found
interface
    subroutine DSTEQR(COMPZ, N, D, E, Z, LDZ, WORK, INFO)
        import :: dp
        character :: COMPZ
        integer :: N, LDZ, INFO
        real(dp) :: D(*), E(*), Z(LDZ, *), WORK(*)
    end subroutine DSTEQR
end interface

contains

end module gjp_lapack

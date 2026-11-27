! BEGIN_HEADER
! -----------------------------------------------------------------------------
! Gauss-Jacobi Quadrature Implementation
! Authors: Rohit Goswami <rgoswami[at]ieee.org>
! Source: GaussJacobiQuad Library
! License: MIT
! GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
! -----------------------------------------------------------------------------
! END_HEADER
!> @brief Module for defining types and precision levels for Gauss-Jacobi Polynomial (GJP) calculations.
!!
!! @details This module defines various real kinds for numerical calculations in single, double, high,
!! and quadruple precisions. It also defines a data type for sparse representation of Jacobi matrices.
!! The module makes these types and precision levels public so they can be used in other modules and programs.
!!
!! Types:
!! - `gjp_sparse_matrix`: A type for representing Jacobi matrices in a sparse format.
!!
!! Precision Levels:
!! - `sp`: Single Precision
!! - `dp`: Double Precision
!! - `hp`: High Precision
!! - `qp`: Quadruple Precision
module gjp_types
use, intrinsic :: iso_c_binding, only: c_int
implicit none
private
public sp, dp, hp, qp, gjp_sparse_matrix, GJP_ERR_BOUNDS

!> @brief Define various kinds for real numbers.
integer, parameter :: dp = kind(0.d0), & ! double precision
                      hp = selected_real_kind(15), & ! high precision
                      qp = selected_real_kind(32), & ! quadruple precision
                      sp = kind(0.) ! single precision
! Shared with GaussJacobiQuad and GaussJacobiQuadCCompat. Matches the C enum.
integer(c_int), parameter :: GJP_ERR_BOUNDS = 7_c_int

!> @brief Sparse representation of a Jacobi matrix.
    !!
    !! @details The type gjp_sparse_matrix represents a Jacobi matrix in a sparse format.
    !! It holds the diagonal and off-diagonal elements of the matrix, as the other elements are zeros.
    !! The sparse format is suitable for efficient storage and computation when dealing with Jacobi matrices.
type gjp_sparse_matrix
    real(dp), allocatable :: diagonal(:) !< Diagonal elements of the Jacobi matrix.
    real(dp), allocatable :: off_diagonal(:) !< Off-diagonal elements of the Jacobi matrix.
end type gjp_sparse_matrix

end module gjp_types

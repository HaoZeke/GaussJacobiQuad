! BEGIN_HEADER
! -----------------------------------------------------------------------------
! Gauss-Jacobi Quadrature Implementation
! Authors: Rohit Goswami <rgoswami[at]ieee.org>
! Source: GaussJacobiQuad Library
! License: MIT
! GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
! Date: 2023-08-26
! Commit: 5dd0ffe
! -----------------------------------------------------------------------------
! This code is part of the GaussJacobiQuad library, providing an efficient
! implementation for Gauss-Jacobi quadrature nodes and weights computation.
! -----------------------------------------------------------------------------
! END_HEADER
!> @brief Compatibility layer for Gauss-Jacobi quadrature with C bindings.
!
!! @details This module provides a C-compatible interface to the Gauss-Jacobi quadrature methods
!! implemented in the GaussJacobiQuad module. It's designed to make the Fortran implementations
!! accessible from C/C++ code. It wraps the original Fortran subroutines with ones that are C-compatible,
!! using the `iso_c_binding` module for the necessary type conversions.
module GaussJacobiQuadCCompat
use iso_c_binding, only: c_double, c_int
use GaussJacobiQuad, only: gauss_jacobi_rec, gauss_jacobi_gw, gauss_jacobi_algo665
implicit none
contains

!> @brief C-compatible wrapper for `gauss_jacobi_rec` subroutine.
!>
!> This subroutine is a C-compatible wrapper that calls the original `gauss_jacobi_rec`
!> subroutine for calculating Gauss-Jacobi quadrature nodes and weights using the recurrence method.
!>
!> @param[in] npts Number of quadrature points.
!> @param[in] alpha Parameter alpha in the Jacobi polynomial. Must be greater than -1.
!> @param[in] beta Parameter beta in the Jacobi polynomial. Must be greater than -1.
!> @param[out] x Quadrature nodes.
!> @param[out] wts Quadrature weights.
subroutine gauss_jacobi_rec_c(npts, alpha, beta, x, wts) bind(C, name="gauss_jacobi_rec_c")
    integer(c_int), intent(in) :: npts
    real(c_double), intent(in) :: alpha, beta
    real(c_double), intent(out) :: x(npts), wts(npts)
    call gauss_jacobi_rec(npts, alpha, beta, x, wts)
end subroutine gauss_jacobi_rec_c

!> @brief C-compatible wrapper for `gauss_jacobi_gw` subroutine.
!>
!> This subroutine is a C-compatible wrapper that calls the original `gauss_jacobi_gw`
!> subroutine for calculating Gauss-Jacobi quadrature nodes and weights using the "gw" method.
!>
!> @param[in] npts Number of quadrature points.
!> @param[in] alpha Parameter alpha in the Jacobi polynomial. Must be greater than -1.
!> @param[in] beta Parameter beta in the Jacobi polynomial. Must be greater than -1.
!> @param[out] x Quadrature nodes.
!> @param[out] wts Quadrature weights.
subroutine gauss_jacobi_gw_c(npts, alpha, beta, x, wts) bind(C, name="gauss_jacobi_gw_c")
    integer(c_int), intent(in) :: npts
    real(c_double), intent(in) :: alpha, beta
    real(c_double), intent(out) :: x(npts), wts(npts)
    call gauss_jacobi_gw(npts, alpha, beta, x, wts)
end subroutine gauss_jacobi_gw_c

!> @brief C-compatible wrapper for `gauss_jacobi_algo665` subroutine.
!>
!> This subroutine is a C-compatible wrapper that calls the original `gauss_jacobi_algo665`
!> subroutine for calculating Gauss-Jacobi quadrature nodes and weights using the "algo665" method.
!>
!> @param[in] npts Number of quadrature points.
!> @param[in] alpha Parameter alpha in the Jacobi polynomial. Must be greater than -1.
!> @param[in] beta Parameter beta in the Jacobi polynomial. Must be greater than -1.
!> @param[out] x Quadrature nodes.
!> @param[out] wts Quadrature weights.
subroutine gauss_jacobi_algo665_c(npts, alpha, beta, x, wts) bind(C, name="gauss_jacobi_algo665_c")
    integer(c_int), intent(in) :: npts
    real(c_double), intent(in) :: alpha, beta
    real(c_double), intent(out) :: x(npts), wts(npts)
    call gauss_jacobi_algo665(npts, alpha, beta, x, wts)
end subroutine gauss_jacobi_algo665_c

end module GaussJacobiQuadCCompat

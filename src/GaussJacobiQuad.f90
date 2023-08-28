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
!> @brief Overall driver for Gauss-Jacobi quadrature.
!
!! @details This module implements the Gauss-Jacobi quadrature method for numerical integration.
!! It provides a subroutine for obtaining the weights and nodes, which dispatches to multiple implementations
!! based on the method provided. Available methods include "rec".
module GaussJacobiQuad
use gjp_rec, only: gauss_jacobi_rec
use gjp_gw, only: gauss_jacobi_gw
use gjp_algo665, only: gauss_jacobi_algo665
use gjp_types, only: dp
implicit none
contains

!> @brief Compute the Gauss-Jacobi quadrature nodes and weights.
!>
!> This subroutine calculates the Gauss-Jacobi quadrature nodes and weights for the given parameters @f$\alpha@f$ and @f$\beta@f$,
!> using the specified method. Gauss-Jacobi quadrature is used to approximate integrals of the form:
!> \[
!>   \int_{-1}^{1} (1 - x)^\alpha (1 + x)^\beta f(x) \,dx \approx \sum_{i=1}^{npts} wts_i f(x_i)
!> \]
!> where the weights and nodes are calculated with the Jacobi polynomial, which is defined as:
!> \[
!>   P_n^{(\alpha, \beta)}(x) = \frac{(\alpha + 1)_n}{n!} \sum_{k=0}^n \binom{n}{k} \frac{(\beta + 1)_{n-k}}{(n-k)!} \left( \frac{x-1}{2} \right)^k \left( \frac{x+1}{2} \right)^{n-k}
!> \]
!>
!> @param[in] npts Number of quadrature points.
!> @param[in] alpha Parameter alpha in the Jacobi polynomial. Must be greater than -1.
!> @param[in] beta Parameter beta in the Jacobi polynomial. Must be greater than -1.
!> @param[out] x Quadrature nodes.
!> @param[out] wts Quadrature weights.
!> @param[in] method Method used for calculation. Supported methods are "rec" and "gw".
subroutine gauss_jacobi(npts, alpha, beta, x, wts, method)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    character(len=:), allocatable, intent(in) :: method

    if (npts <= 0) then
        error stop "Number of points must be positive"
    end if

    if (alpha <= -1.0_dp) then
        error stop "alpha must be greater than -1"
    end if

    if (beta <= -1.0_dp) then

        error stop "beta must be greater than -1"
    end if

    select case (trim(method))
    case ("rec") ! Fails at high beta
        call gauss_jacobi_rec(npts, alpha, beta, x, wts)
    case ("gw") ! Accurate for high beta
        call gauss_jacobi_gw(npts, alpha, beta, x, wts)
    case ("algo665")
        call gauss_jacobi_algo665(npts, alpha, beta, x, wts)
    case default
        print*,"Error: Unknown method specified:", method
        print*,"Supported methods: 'rec', 'gw', 'algo665'"
        error stop
    end select
end subroutine gauss_jacobi

end module

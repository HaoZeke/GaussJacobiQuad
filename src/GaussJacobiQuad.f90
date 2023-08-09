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
!> @brief Overall driver for Gauss-Jacobi quadrature.
!
!! @details This module implements the Gauss-Jacobi quadrature method for numerical integration.
!! It provides a subroutine for obtaining the weights and nodes, which dispatches to multiple implementations
!! based on the method provided. Available methods include "recurrence".
module GaussJacobiQuad
use gjp_rec, only: gauss_jacobi_rec
use gjp_gw, only: gauss_jacobi_gw
use gjp_types, only: dp
implicit none
contains

subroutine gauss_jacobi(npts, alpha, beta, x, wts, method)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    character(len=:), allocatable, intent(in) :: method

    if (alpha <= -1.0_dp) then
        print*,"Error: alpha must be greater than -1"
        error stop
    end if

    if (beta <= -1.0_dp) then
        print*,"Error: beta must be greater than -1"
        error stop
    end if

    select case (trim(method))
    case ("recurrence") ! Fails at high beta
        call gauss_jacobi_rec(npts, alpha, beta, x, wts)
    case ("gw") ! Accurate for high beta
        call gauss_jacobi_gw(npts, alpha, beta, x, wts)
    case default
        print*,"Error: Unknown method specified:", method
        print*,"Supported methods: 'recurrence', 'gw''"
        error stop
    end select
end subroutine gauss_jacobi

end module

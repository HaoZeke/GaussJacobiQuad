! BEGIN_HEADER
! -----------------------------------------------------------------------------
! Gauss-Jacobi Quadrature Implementation
! Authors: Rohit Goswami <rgoswami[at]ieee.org>
! Source: GaussJacobiQuad Library
! License: MIT
! GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
! Date: 2023-08-25
! Commit: f68c09f
! -----------------------------------------------------------------------------
! This code is part of the GaussJacobiQuad library, providing an efficient
! implementation for Gauss-Jacobi quadrature nodes and weights computation.
! -----------------------------------------------------------------------------
! END_HEADER
module GaussJacobiQuadCCompat
use iso_c_binding, only: c_double, c_int
use GaussJacobiQuad, only: gauss_jacobi_rec, gauss_jacobi_gw
implicit none
contains

subroutine gauss_jacobi_rec_c(npts, alpha, beta, x, wts) bind(C, name="gauss_jacobi_rec_c")
    integer(c_int), intent(in) :: npts
    real(c_double), intent(in) :: alpha, beta
    real(c_double), intent(out) :: x(npts), wts(npts)
    call gauss_jacobi_rec(npts, alpha, beta, x, wts)
end subroutine gauss_jacobi_rec_c

subroutine gauss_jacobi_gw_c(npts, alpha, beta, x, wts) bind(C, name="gauss_jacobi_gw_c")
    integer(c_int), intent(in) :: npts
    real(c_double), intent(in) :: alpha, beta
    real(c_double), intent(out) :: x(npts), wts(npts)
    call gauss_jacobi_gw(npts, alpha, beta, x, wts)
end subroutine gauss_jacobi_gw_c

end module GaussJacobiQuadCCompat

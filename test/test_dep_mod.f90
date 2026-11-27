! BEGIN_HEADER
! -----------------------------------------------------------------------------
! Gauss-Jacobi Quadrature Implementation
! Authors: Rohit Goswami <rgoswami[at]ieee.org>
! Source: GaussJacobiQuad Library
! License: MIT
! GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
! -----------------------------------------------------------------------------
! END_HEADER
program test_dep_mod
use GaussJacobiQuad, only: dp
implicit none
real(dp) :: x
x = 1.0_dp
if (x <= 0.0_dp) error stop 1
end program test_dep_mod

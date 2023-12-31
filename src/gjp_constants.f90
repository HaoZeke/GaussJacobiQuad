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
!> @brief Constants contain more digits than double precision, so that
!> they are rounded correctly.
!> @details Single letter constants contain underscore so
!> that they do not clash with user variables ("e" and "i" are frequently used as
!> loop variables). The public parameters are enumerated below:
!> @param c_year values correspond to different CODATA standards for the inverse fine-structure constant
!> @param Ha2eV_year values is the conversion factor for converting Hartree to eV as per varying CODATA standards
!> @param e_ is Euler's constant
!> @param i_ denotes one unit length along the complex axis as a conversion helper

module gjp_constants
use gjp_types, only: dp
implicit none
private
public pi, e_, i_, Ha2eV_2010, c_2010, c_2006, c_1998, c_1986

real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp
real(dp), parameter :: e_ = 2.7182818284590452353602874713527_dp
complex(dp), parameter :: i_ = (0, 1)

real(dp), parameter :: Ha2eV_2010 = 27.21138505_dp ! 1 Ha = (1 * Ha2eV) eV
! Standard uncertainty:      0.000000000053 eV (Source: 2010 CODATA)

! speed of light in atomic units  --> inverse fine-structure constant
real(dp), parameter :: c_2010 = 137.035999074_dp
real(dp), parameter :: c_2006 = 137.035999679_dp
real(dp), parameter :: c_2002 = 137.03599911_dp
real(dp), parameter :: c_1998 = 137.03599976_dp
real(dp), parameter :: c_1986 = 137.0359895_dp ! To compare with dftatom
end module gjp_constants

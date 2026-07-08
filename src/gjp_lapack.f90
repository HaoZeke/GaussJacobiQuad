! BEGIN_HEADER
! -----------------------------------------------------------------------------
! Gauss-Jacobi Quadrature Implementation
! Authors: Rohit Goswami <rgoswami[at]ieee.org>
! Source: GaussJacobiQuad Library
! License: MIT
! GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
! -----------------------------------------------------------------------------
! END_HEADER
!> @brief LAPACK Interface for Gauss Jacobi Polynomials
!> @details Portable wrapper around the tridiagonal eigenproblem via LAPACKE
!> (C ABI) to avoid gfortran/OpenBLAS Fortran character-string ABI segfaults.
module gjp_lapack
use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char
implicit none

integer, parameter :: dp = kind(0.d0)
integer(c_int), parameter, private :: LAPACK_COL_MAJOR = 102

interface
    function lapacke_dsteqr(matrix_layout, compz, n, d, e, z, ldz) &
        bind(C, name="LAPACKE_dsteqr")
        import :: c_int, c_double, c_char
        integer(c_int) :: lapacke_dsteqr
        integer(c_int), value :: matrix_layout
        character(kind=c_char), value :: compz
        integer(c_int), value :: n, ldz
        real(c_double) :: d(*), e(*), z(*)
    end function lapacke_dsteqr
end interface

contains

!> Fortran-shaped wrapper matching the historical DSTEQR call site in gjp_gw.
subroutine DSTEQR(COMPZ, N, D, E, Z, LDZ, WORK, INFO)
    character(len=1), intent(in) :: COMPZ
    integer, intent(in) :: N, LDZ
    integer, intent(out) :: INFO
    real(dp), intent(inout) :: D(*), E(*), Z(LDZ, *), WORK(*)
    integer(c_int) :: n_c, ldz_c, info_c
    character(kind=c_char) :: compz_c

    ! WORK is unused by LAPACKE (allocates internally)
    associate (unused => WORK(1))
    end associate

    n_c = int(N, c_int)
    ldz_c = int(LDZ, c_int)
    compz_c = COMPZ
    info_c = lapacke_dsteqr(LAPACK_COL_MAJOR, compz_c, n_c, D, E, Z, ldz_c)
    INFO = int(info_c)
end subroutine DSTEQR

end module gjp_lapack

! BEGIN_HEADER
! -----------------------------------------------------------------------------
! Gauss-Jacobi Quadrature Implementation
! Authors: Rohit Goswami <rgoswami[at]ieee.org>
! Source: GaussJacobiQuad Library
! License: MIT
! GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
! -----------------------------------------------------------------------------
! END_HEADER
!> @brief Overall driver for Gauss-Jacobi quadrature with multi-algorithm dispatch.
module GaussJacobiQuad
use gjp_rec, only: gauss_jacobi_rec, gauss_jacobi_rec_caf
use gjp_gw, only: gauss_jacobi_gw
use gjp_algo665, only: gauss_jacobi_algo665, gauss_jacobi_algo665_dc
use gjp_sturm, only: gauss_jacobi_sturm, gauss_jacobi_sturm_caf
use gjp_bogaert, only: gauss_jacobi_bogaert, gauss_jacobi_bogaert_caf
use gjp_glr, only: gauss_jacobi_glr, gauss_jacobi_glr_caf
use gjp_caf, only: gauss_jacobi_batch_caf
use gjp_auto, only: select_method_auto
use gjp_types, only: dp
implicit none
public :: gauss_jacobi, gauss_jacobi_rule, gauss_jacobi_batch_caf
public :: select_method_auto
contains

!> Legacy/string dispatch. Prefer =gauss_jacobi_rule= for optional method/auto.
subroutine gauss_jacobi(npts, alpha, beta, x, wts, method)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    character(len=:), allocatable, intent(in) :: method
    call gauss_jacobi_rule(npts, alpha, beta, x, wts, method)
end subroutine gauss_jacobi

!> Single public entry: optional method; blank / "auto" selects a policy.
!>
!> @param method Optional. If absent, empty, or "auto", uses select_method_auto
!>        (CAF names when num_images()>1). Forced names include *_caf variants.
subroutine gauss_jacobi_rule(npts, alpha, beta, x, wts, method)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    character(len=*), intent(in), optional :: method
    character(len=:), allocatable :: m

    if (npts <= 0) error stop "Number of points must be positive"
    if (alpha <= -1.0_dp) error stop "alpha must be greater than -1"
    if (beta <= -1.0_dp) error stop "beta must be greater than -1"

    if (present(method)) then
        if (len_trim(method) == 0 .or. trim(method) == "auto") then
            m = select_method_auto(npts, alpha, beta)
        else
            m = trim(method)
        end if
    else
        m = select_method_auto(npts, alpha, beta)
    end if

    select case (m)
    case ("rec")
        call gauss_jacobi_rec(npts, alpha, beta, x, wts)
    case ("rec_caf")
        call gauss_jacobi_rec_caf(npts, alpha, beta, x, wts)
    case ("gw", "gw_caf")
        call gauss_jacobi_gw(npts, alpha, beta, x, wts)
    case ("algo665", "algo665_caf")
        call gauss_jacobi_algo665(npts, alpha, beta, x, wts)
    case ("algo665_dc")
        call gauss_jacobi_algo665_dc(npts, alpha, beta, x, wts, use_caf=.false.)
    case ("algo665_dc_caf")
        call gauss_jacobi_algo665_dc(npts, alpha, beta, x, wts, use_caf=.true.)
    case ("sturm")
        call gauss_jacobi_sturm(npts, alpha, beta, x, wts)
    case ("sturm_caf")
        call gauss_jacobi_sturm_caf(npts, alpha, beta, x, wts)
    case ("glr")
        call gauss_jacobi_glr(npts, alpha, beta, x, wts)
    case ("glr_caf")
        call gauss_jacobi_glr_caf(npts, alpha, beta, x, wts)
    case ("bogaert")
        call gauss_jacobi_bogaert(npts, alpha, beta, x, wts)
    case ("bogaert_caf")
        call gauss_jacobi_bogaert_caf(npts, alpha, beta, x, wts)
    case default
        print*,"Error: Unknown method:", m
        print*,"Supported: auto, rec, rec_caf, gw, algo665, algo665_dc, algo665_dc_caf,"
        print*,"           sturm, sturm_caf, glr, glr_caf, bogaert, bogaert_caf"
        error stop
    end select
end subroutine gauss_jacobi_rule

end module GaussJacobiQuad

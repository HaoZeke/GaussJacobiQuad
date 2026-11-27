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
use iso_c_binding, only: c_double, c_int
implicit none
public :: gauss_jacobi, gauss_jacobi_rule, gauss_jacobi_batch_caf
public :: select_method_auto, gauss_jacobi_scale
private :: scale_interval, gjp_scale_f
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

!> Map Gauss-Jacobi nodes and weights from [-1, 1] onto [a, b].
!>
!> Nodes go to t = (b-a)/2 * x + (a+b)/2.
!> For the Jacobi weight (1-x)**alpha * (1+x)**beta, weights are multiplied
!> by ((b-a)/2)**(alpha+beta+1) so the mapped rule integrates
!> (b-t)**alpha * (t-a)**beta * f(t) on [a, b].
!> Requires finite a < b and finite alpha, beta > -1.
!> (a+b) and (b-a) are not formed when that would overflow.
!> A non-finite scaled node or weight is an error and does not write.
subroutine gauss_jacobi_scale(npts, alpha, beta, a, b, x, wts)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta, a, b
    real(dp), intent(inout) :: x(npts), wts(npts)
    integer :: st

    st = scale_interval(npts, alpha, beta, a, b, x, wts)
    if (st == 1) error stop "Number of points must be positive"
    if (st == 2) error stop "alpha must be greater than -1"
    if (st == 3) error stop "beta must be greater than -1"
    if (st == 7) error stop "invalid interval or scaled weights overflow"
end subroutine gauss_jacobi_scale

!> C ABI. Status: 0 ok, 1 npts, 2 alpha, 3 beta, 7 bounds (GJP_ERR_BOUNDS).
function gjp_scale_f(npts, alpha, beta, a, b, x, wts) result(status) &
    bind(C, name="gjp_scale_f")
    integer(c_int), intent(in), value :: npts
    real(c_double), intent(in), value :: alpha, beta, a, b
    real(c_double), intent(inout) :: x(npts), wts(npts)
    integer(c_int) :: status
    real(dp), allocatable :: xd(:), wd(:)
    integer :: st

    if (npts <= 0_c_int) then
        status = 1_c_int
        return
    end if
    allocate (xd(npts), wd(npts))
    xd = real(x, dp)
    wd = real(wts, dp)
    st = scale_interval(int(npts), real(alpha, dp), real(beta, dp), &
                        real(a, dp), real(b, dp), xd, wd)
    status = int(st, c_int)
    if (status == 0_c_int) then
        x = real(xd, c_double)
        wts = real(wd, c_double)
    end if
end function gjp_scale_f

integer function scale_interval(npts, alpha, beta, a, b, x, wts) result(status)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta, a, b
    real(dp), intent(inout) :: x(npts), wts(npts)
    real(dp) :: half, fac, mid, t, tw
    real(dp), allocatable :: xt(:), wt(:)
    integer :: i
    logical :: lease_ok, repair

    status = 0
    if (npts <= 0) then
        status = 1
        return
    end if
    if (.not. ieee_is_finite(alpha) .or. alpha <= -1.0_dp) then
        status = 2
        return
    end if
    if (.not. ieee_is_finite(beta) .or. beta <= -1.0_dp) then
        status = 3
        return
    end if
    if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b) .or. .not. (a < b)) then
        status = 7
        return
    end if
    ! Lease map when 0.5*(b-a) and 0.5*(a+b) are finite. A node whose
    ! image is non-finite, outside [a, b], or collapsed onto an endpoint
    ! uses (a/2)*(1-x)+(b/2)*(1+x) and is clamped into [a, b].
    half = 0.5_dp * (b - a)
    mid = 0.5_dp * (a + b)
    lease_ok = ieee_is_finite(half) .and. ieee_is_finite(mid)
    if (.not. lease_ok) then
        half = 0.5_dp * b - 0.5_dp * a
        mid = 0.5_dp * a + 0.5_dp * b
        if (.not. ieee_is_finite(half)) then
            status = 7
            return
        end if
    end if
    fac = half**(alpha + beta + 1.0_dp)
    if (.not. ieee_is_finite(fac)) then
        status = 7
        return
    end if
    allocate (xt(npts), wt(npts))
    do i = 1, npts
        if (lease_ok) then
            t = half * x(i) + mid
        else
            t = mid + half * x(i)
        end if
        repair = .not. (ieee_is_finite(t) .and. t >= a .and. t <= b)
        if (ieee_is_finite(t) .and. t == a .and. x(i) /= -1.0_dp) repair = .true.
        if (ieee_is_finite(t) .and. t == b .and. x(i) /= 1.0_dp) repair = .true.
        if (repair) then
            t = (0.5_dp * a) * (1.0_dp - x(i)) + (0.5_dp * b) * (1.0_dp + x(i))
            if (.not. ieee_is_finite(t)) then
                status = 7
                return
            end if
            if (t < a) t = a
            if (t > b) t = b
        end if
        tw = fac * wts(i)
        if (.not. ieee_is_finite(tw)) then
            status = 7
            return
        end if
        xt(i) = t
        wt(i) = tw
    end do
    x = xt
    wts = wt
end function scale_interval

end module GaussJacobiQuad

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
use gjp_types, only: dp, GJP_ERR_BOUNDS
use iso_c_binding, only: c_double, c_int
use, intrinsic :: iso_fortran_env, only: int64
! Apps link this facade. fpm drops a src module that no executable
! dependency uses, so install would omit gjp_version.mod and its object.
use gjp_version, only: gjp_version_string
implicit none
private :: gjp_version_string
public :: gauss_jacobi, gauss_jacobi_rule, gauss_jacobi_batch_caf
public :: select_method_auto, gauss_jacobi_scale
private :: scale_interval, gjp_scale_f, affine_image
private :: aff_from_double, aff_uadd, aff_usub, aff_ucmp, aff_sadd, aff_smul
private :: aff_bit, aff_sticky, aff_round, pow2_exp
private :: split_dp, z_set_int, z_trim, z_mul_int, z_mul_big, z_pow, z_bitlen, z_bit
private :: z_sticky_below, emit, z_round, round_int, fp_zero, fp_set_one, fp_zero_p
private :: fp_low_nonzero, fp_add_small, fp_add, fp_shl1, fp_mul_raw, fp_take, fp_div_small
private :: fp_ln2, fp_exp, fp_as_z, fp_ge, fp_sub, fp_recip_unit, fp_div_i64, fp_mul_i64
private :: fp_div_i64_small, z_bitlen_small, place_i64_at, split_scaled, read_shifted
private :: frac_from_shift, scale_clamp_mant, factor_hpow, apply_sig, int_pow_fits
private :: factor_int, apply_int, scaled_weight
! Each operand is an integer significand times 2**1074. aff_round divides
! that product by 2**2149. 280 base-2**16 limbs cover every finite binary64
! product (the leading limb stays near 262).
integer, parameter, private :: aff_nl = 280
integer, parameter, private :: aff_lw = 16
integer(int64), parameter, private :: aff_mask = 65535_int64
! Limb buffer for an exact integer power. 5**3181 fits; a 53-bit odd
! base fits through p about 1200. Larger integer powers use the same
! fixed-point log/exp as a non-integer power, rounded once.
integer, parameter, private :: zmax = 4096
integer, parameter, private :: fb = 1280
integer, parameter, private :: fl = 96
integer(int64), parameter, private :: zmask = 65535_int64
character(len=*), parameter, private :: gjp_version_keep = gjp_version_string
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
!> Nodes go to t = (b-a)/2 * x + (a+b)/2, within 1 ulp of that value,
!> and exactly a or b at x = -1 or x = 1.
!> For the Jacobi weight (1-x)**alpha * (1+x)**beta, weights are multiplied
!> by ((b-a)/2)**(alpha+beta+1) so the mapped rule integrates
!> (b-t)**alpha * (t-a)**beta * f(t) on [a, b].
!> Requires finite a < b and finite alpha, beta > -1.
!> (a+b) and (b-a) are not formed when that would overflow.
!> A non-finite scaled node or weight is an error and does not write.
!> An overflowing weight factor is not an error when the product is finite.
subroutine gauss_jacobi_scale(npts, alpha, beta, a, b, x, wts)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta, a, b
    real(dp), intent(inout) :: x(npts), wts(npts)
    integer :: st

    st = scale_interval(npts, alpha, beta, a, b, x, wts)
    if (st == 1) error stop "Number of points must be positive"
    if (st == 2) error stop "alpha must be greater than -1"
    if (st == 3) error stop "beta must be greater than -1"
    if (st == GJP_ERR_BOUNDS) error stop "invalid interval or scaled weights overflow"
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
    real(dp) :: half, mid, p, t, tw, fac
    real(dp), allocatable :: xt(:), wt(:)
    integer :: i, k, wmode, iscal, sigst, nmag
    integer(int64) :: pk, signf, fexp
    integer(int64) :: wsig(fl)
    integer(int64), allocatable :: pmag(:)
    logical :: lease_ok

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
        status = GJP_ERR_BOUNDS
        return
    end if
    ! Weight factor (b-a)/2. Split it when a+b or b-a overflows.
    ! Nodes are affine_image, not this half/mid form.
    ! 0.5*b - 0.5*a is finite for every finite a < b.
    half = 0.5_dp * (b - a)
    mid = 0.5_dp * (a + b)
    lease_ok = ieee_is_finite(half) .and. ieee_is_finite(mid)
    if (.not. lease_ok) half = 0.5_dp * b - 0.5_dp * a
    p = alpha + beta + 1.0_dp
    if (half < 0.0_dp) then
        status = GJP_ERR_BOUNDS
        return
    end if
    ! -ffast-math and -ffinite-math-only are unsupported: they delete
    ! ieee_is_finite, so a non-finite node or weight is not rejected.
    ! Weight plan is independent of the node. p = 1 is one multiply.
    ! A power-of-two half and an integer p stay on scale(), clamped
    ! to ±4096, which matches ldexp across the binary64 range.
    wmode = 0
    iscal = 0
    sigst = 0
    nmag = 0
    signf = 0_int64
    fexp = 0_int64
    fac = 1.0_dp
    allocate (pmag(zmax))
    if (.not. ieee_is_finite(p)) then
        wmode = 3
    else if (p /= 0.0_dp .and. half /= 1.0_dp) then
        k = pow2_exp(half)
        if (k /= -100000 .and. p == aint(p) .and. abs(p) <= 2.0_dp**52) then
            pk = int(p, int64) * int(k, int64)
            if (pk > 4096_int64) then
                iscal = 4096
            else if (pk < -4096_int64) then
                iscal = -4096
            else
                iscal = int(pk)
            end if
            wmode = 1
        else if (p == 1.0_dp) then
            wmode = 2
            fac = half
        else if (p == aint(p) .and. p > 0.0_dp .and. p <= 2.0_dp**52 &
                 .and. int_pow_fits(half, p)) then
            wmode = 4
            call factor_int(int(p, int64), half, fexp, pmag, nmag)
        else
            wmode = 5
            call factor_hpow(p, half, sigst, signf, wsig)
        end if
    end if
    allocate (xt(npts), wt(npts))
    do i = 1, npts
        if (x(i) == -1.0_dp) then
            t = a
        else if (x(i) == 1.0_dp) then
            t = b
        else
            t = affine_image(a, b, x(i))
            if (.not. ieee_is_finite(t)) then
                status = GJP_ERR_BOUNDS
                return
            end if
        end if
        if (.not. ieee_is_finite(t) .or. .not. ieee_is_finite(wts(i))) then
            status = GJP_ERR_BOUNDS
            return
        end if
        call scaled_weight(wts(i), p, half, wmode, iscal, fac, sigst, signf, &
                           wsig, fexp, pmag, nmag, tw)
        if (.not. ieee_is_finite(tw)) then
            status = GJP_ERR_BOUNDS
            return
        end if
        xt(i) = t
        wt(i) = tw
    end do
    x = xt
    wts = wt
end function scale_interval

!> Correctly rounded a + (x+1)*(b-a)/2 (round ties to even).
pure function affine_image(a, b, x) result(t)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(dp), intent(in) :: a, b, x
    real(dp) :: t
    integer(int64) :: ma(0:aff_nl - 1), mb(0:aff_nl - 1), mx(0:aff_nl - 1)
    integer(int64) :: one(0:aff_nl - 1), m1(0:aff_nl - 1), p1(0:aff_nl - 1)
    integer(int64) :: t1(0:aff_nl - 1), t2(0:aff_nl - 1), ms(0:aff_nl - 1)
    logical :: sa, sb, sx, s1, sp, st1, st2, ss

    if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b) .or. .not. ieee_is_finite(x)) then
        t = a + b + x
        return
    end if
    call aff_from_double(a, sa, ma)
    call aff_from_double(b, sb, mb)
    call aff_from_double(x, sx, mx)
    one = 0_int64
    ! 2**1074 sits in limb 67, bit 2 (1074 = 67*16 + 2).
    one(67) = 4_int64
    call aff_sadd(.false., one,.not. sx, mx, s1, m1)
    call aff_sadd(.false., one, sx, mx, sp, p1)
    call aff_smul(sa, ma, s1, m1, st1, t1)
    call aff_smul(sb, mb, sp, p1, st2, t2)
    call aff_sadd(st1, t1, st2, t2, ss, ms)
    t = aff_round(ss, ms)
end function affine_image

pure subroutine aff_from_double(val, neg, mag)
    real(dp), intent(in) :: val
    logical, intent(out) :: neg
    integer(int64), intent(out) :: mag(0:)
    integer(int64) :: bits, frac, sig, rest, chunk, cur, carry
    integer :: expf, shift, pos, limb, off, i

    mag = 0_int64
    neg = .false.
    if (iand(transfer(val, 0_int64), not(ishft(1_int64, 63))) == 0_int64) return
    bits = transfer(val, 0_int64)
    neg = btest(bits, 63)
    expf = int(iand(ishft(bits, -52), 2047_int64))
    frac = iand(bits, 4503599627370495_int64)
    if (expf == 0) then
        sig = frac
        shift = 0
    else
        sig = frac + 4503599627370496_int64
        shift = expf - 1
    end if
    rest = sig
    pos = shift
    do while (rest /= 0_int64)
        chunk = iand(rest, aff_mask)
        rest = ishft(rest, -16)
        limb = pos / aff_lw
        off = mod(pos, aff_lw)
        cur = mag(limb) + ishft(chunk, off)
        mag(limb) = iand(cur, aff_mask)
        carry = ishft(cur, -16)
        i = limb + 1
        do while (carry /= 0_int64)
            cur = mag(i) + carry
            mag(i) = iand(cur, aff_mask)
            carry = ishft(cur, -16)
            i = i + 1
        end do
        pos = pos + aff_lw
    end do
end subroutine aff_from_double

pure subroutine aff_uadd(aa, ba, c)
    integer(int64), intent(in) :: aa(0:), ba(0:)
    integer(int64), intent(out) :: c(0:)
    integer(int64) :: carry, cur
    integer :: i

    carry = 0_int64
    do i = 0, size(aa) - 1
        cur = aa(i) + ba(i) + carry
        c(i) = iand(cur, aff_mask)
        carry = ishft(cur, -16)
    end do
end subroutine aff_uadd

pure subroutine aff_usub(aa, ba, c)
    integer(int64), intent(in) :: aa(0:), ba(0:)
    integer(int64), intent(out) :: c(0:)
    integer(int64) :: cur, borrow
    integer :: i

    borrow = 0_int64
    do i = 0, size(aa) - 1
        cur = aa(i) - ba(i) - borrow
        if (cur < 0_int64) then
            cur = cur + 65536_int64
            borrow = 1_int64
        else
            borrow = 0_int64
        end if
        c(i) = cur
    end do
end subroutine aff_usub

pure function aff_ucmp(aa, ba) result(cmp)
    integer(int64), intent(in) :: aa(0:), ba(0:)
    integer :: cmp
    integer :: i

    cmp = 0
    do i = size(aa) - 1, 0, -1
        if (aa(i) /= ba(i)) then
            if (aa(i) > ba(i)) then
                cmp = 1
            else
                cmp = -1
            end if
            return
        end if
    end do
end function aff_ucmp

pure subroutine aff_sadd(sa, ma, sb, mb, so, mo)
    logical, intent(in) :: sa, sb
    integer(int64), intent(in) :: ma(0:), mb(0:)
    logical, intent(out) :: so
    integer(int64), intent(out) :: mo(0:)
    integer :: cmp

    if (sa .eqv. sb) then
        so = sa
        call aff_uadd(ma, mb, mo)
        if (all(mo == 0_int64)) so = .false.
    else
        cmp = aff_ucmp(ma, mb)
        if (cmp == 0) then
            so = .false.
            mo = 0_int64
        else if (cmp > 0) then
            so = sa
            call aff_usub(ma, mb, mo)
        else
            so = sb
            call aff_usub(mb, ma, mo)
        end if
    end if
end subroutine aff_sadd

pure subroutine aff_smul(sa, ma, sb, mb, so, mo)
    logical, intent(in) :: sa, sb
    integer(int64), intent(in) :: ma(0:), mb(0:)
    logical, intent(out) :: so
    integer(int64), intent(out) :: mo(0:)
    integer(int64) :: acc, prod
    integer :: i, j, n

    if (all(ma == 0_int64) .or. all(mb == 0_int64)) then
        so = .false.
        mo = 0_int64
        return
    end if
    so = sa .neqv. sb
    n = size(ma)
    mo = 0_int64
    do i = 0, n - 1
        if (ma(i) == 0_int64) cycle
        acc = 0_int64
        do j = 0, n - 1 - i
            prod = ma(i) * mb(j) + mo(i + j) + acc
            mo(i + j) = iand(prod, aff_mask)
            acc = ishft(prod, -16)
        end do
    end do
end subroutine aff_smul

pure function aff_bit(mag, pos) result(on)
    integer(int64), intent(in) :: mag(0:)
    integer, intent(in) :: pos
    logical :: on
    integer :: limb, off

    if (pos < 0) then
        on = .false.
        return
    end if
    limb = pos / aff_lw
    off = mod(pos, aff_lw)
    on = btest(mag(limb), off)
end function aff_bit

pure function aff_sticky(mag, pos) result(anyb)
    integer(int64), intent(in) :: mag(0:)
    integer, intent(in) :: pos
    logical :: anyb
    integer :: limb, off, i
    integer(int64) :: low

    anyb = .false.
    if (pos <= 0) return
    limb = pos / aff_lw
    off = mod(pos, aff_lw)
    do i = 0, limb - 1
        if (mag(i) /= 0_int64) then
            anyb = .true.
            return
        end if
    end do
    if (off > 0) then
        low = ishft(1_int64, off) - 1_int64
        if (iand(mag(limb), low) /= 0_int64) anyb = .true.
    end if
end function aff_sticky

pure function aff_round(neg, mag) result(y)
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
    logical, intent(in) :: neg
    integer(int64), intent(in) :: mag(0:)
    real(dp) :: y
    ! mag is the image times 2**2149 (2**1074 from each factor).
    integer, parameter :: sh = 2149
    integer, parameter :: sub = 1075
    integer :: i, k, e, j
    integer(int64) :: sig, bits, quot
    logical :: rb, sticky, lsb

    k = -1
    do i = size(mag) - 1, 0, -1
        if (mag(i) /= 0_int64) then
            k = i * aff_lw + (63 - leadz(mag(i)))
            exit
        end if
    end do
    if (k < 0) then
        y = 0.0_dp
        return
    end if
    e = k - sh
    if (e >= 1024) then
        y = ieee_value(0.0_dp, ieee_positive_inf)
        if (neg) y = -y
        return
    end if
    if (e < -1022) then
        quot = 0_int64
        if (k >= sub) then
            do j = sub, k
                if (aff_bit(mag, j)) quot = ior(quot, ishft(1_int64, j - sub))
            end do
        end if
        rb = aff_bit(mag, sub - 1)
        sticky = aff_sticky(mag, sub - 1)
        lsb = btest(quot, 0)
        if (rb .and. (sticky .or. lsb)) quot = quot + 1_int64
        bits = quot
        if (neg) bits = ior(bits, ishft(1_int64, 63))
        y = transfer(bits, 0.0_dp)
        return
    end if
    sig = 0_int64
    do j = 0, 52
        if (aff_bit(mag, k - j)) sig = ior(sig, ishft(1_int64, 52 - j))
    end do
    rb = aff_bit(mag, k - 53)
    sticky = aff_sticky(mag, k - 53)
    lsb = btest(sig, 0)
    if (rb .and. (sticky .or. lsb)) sig = sig + 1_int64
    if (sig == ishft(1_int64, 53)) then
        sig = ishft(sig, -1)
        e = e + 1
    end if
    if (e >= 1024) then
        y = ieee_value(0.0_dp, ieee_positive_inf)
        if (neg) y = -y
        return
    end if
    bits = ior(ishft(int(e + 1023, int64), 52), sig - ishft(1_int64, 52))
    if (neg) bits = ior(bits, ishft(1_int64, 63))
    y = transfer(bits, 0.0_dp)
end function aff_round

!> Exponent k with h = 2**k, or -100000 when h is not a positive power of two.
!> h = 0 stands for the flushed half-width 2**-1075.
pure function pow2_exp(h) result(k)
    real(dp), intent(in) :: h
    integer :: k
    integer(int64) :: bits, frac
    integer :: expf, b
    k = -100000
    if (h == 0.0_dp) then
        k = -1075
        return
    end if
    if (h < 0.0_dp) return
    bits = transfer(h, 0_int64)
    expf = int(iand(ishft(bits, -52), 2047_int64))
    frac = iand(bits, 4503599627370495_int64)
    if (expf == 0) then
        if (frac /= 0_int64 .and. iand(frac, frac - 1_int64) == 0_int64) then
            b = trailz(frac)
            k = b - 1074
        end if
    else if (frac == 0_int64) then
        k = expf - 1023
    end if
end function pow2_exp

subroutine split_dp(x, mant, exp)
    real(dp), intent(in) :: x
    integer(int64), intent(out) :: mant, exp
    integer(int64) :: bits, frac
    integer :: expf
    mant = 0_int64
    exp = 0_int64
    if (x == 0.0_dp) return
    bits = transfer(abs(x), 0_int64)
    expf = int(iand(ishft(bits, -52), 2047_int64))
    frac = iand(bits, 4503599627370495_int64)
    if (expf == 0) then
        mant = frac
        exp = -1074_int64
    else
        mant = ior(frac, ishft(1_int64, 52))
        exp = int(expf, int64) - 1023_int64 - 52_int64
    end if
end subroutine split_dp

subroutine z_set_int(d, n, v)
    integer(int64), intent(out) :: d(:)
    integer, intent(out) :: n
    integer(int64), intent(in) :: v
    integer(int64) :: x
    n = 0
    x = v
    do while (x > 0_int64)
        n = n + 1
        d(n) = iand(x, zmask)
        x = ishft(x, -16)
    end do
end subroutine z_set_int

subroutine z_trim(d, n)
    integer(int64), intent(in) :: d(:)
    integer, intent(inout) :: n
    do while (n > 0)
        if (d(n) /= 0_int64) exit
        n = n - 1
    end do
end subroutine z_trim

subroutine z_mul_int(d, n, m)
    integer(int64), intent(inout) :: d(:)
    integer, intent(inout) :: n
    integer(int64), intent(in) :: m
    integer(int64) :: parts(4), acc, piece, tmp(zmax)
    integer :: s, i, ntmp, k
    if (m == 0_int64 .or. n == 0) then
        n = 0
        return
    end if
    if (m == 1_int64) return
    parts(1) = iand(m, zmask)
    parts(2) = iand(ishft(m, -16), zmask)
    parts(3) = iand(ishft(m, -32), zmask)
    parts(4) = iand(ishft(m, -48), zmask)
    tmp = 0_int64
    ntmp = 0
    do s = 1, 4
        if (parts(s) == 0_int64) cycle
        acc = 0_int64
        do i = 1, n
            piece = d(i) * parts(s) + acc + tmp(i + s - 1)
            tmp(i + s - 1) = iand(piece, zmask)
            acc = ishft(piece, -16)
        end do
        k = n + s
        do while (acc > 0_int64)
            piece = tmp(k) + acc
            tmp(k) = iand(piece, zmask)
            acc = ishft(piece, -16)
            k = k + 1
        end do
        if (k - 1 > ntmp) ntmp = k - 1
    end do
    n = ntmp
    do i = 1, n
        d(i) = tmp(i)
    end do
end subroutine z_mul_int

subroutine z_mul_big(a, na, b, nb, c, nc)
    integer(int64), intent(in) :: a(:), b(:)
    integer, intent(in) :: na, nb
    integer(int64), intent(out) :: c(:)
    integer, intent(out) :: nc
    integer(int64) :: acc, piece
    integer :: i, j, k
    c = 0_int64
    nc = 0
    if (na == 0 .or. nb == 0) return
    do i = 1, na
        acc = 0_int64
        do j = 1, nb
            piece = c(i + j - 1) + a(i) * b(j) + acc
            c(i + j - 1) = iand(piece, zmask)
            acc = ishft(piece, -16)
        end do
        k = i + nb
        do while (acc > 0_int64)
            if (k > size(c)) error stop "z_mul_big overflow"
            piece = c(k) + acc
            c(k) = iand(piece, zmask)
            acc = ishft(piece, -16)
            k = k + 1
        end do
    end do
    nc = na + nb + 2
    call z_trim(c, nc)
end subroutine z_mul_big

subroutine z_pow(base, p, out, nout)
    integer(int64), intent(in) :: base, p
    integer(int64), intent(out) :: out(:)
    integer, intent(out) :: nout
    integer(int64), allocatable :: a(:), r(:), t(:)
    integer :: na, nr, nt
    integer(int64) :: pp
    allocate (a(zmax), r(zmax), t(zmax))
    call z_set_int(a, na, base)
    call z_set_int(r, nr, 1_int64)
    pp = p
    do while (pp > 0_int64)
        if (iand(pp, 1_int64) == 1_int64) then
            call z_mul_big(r, nr, a, na, t, nt)
            nr = nt
            r(1:nr) = t(1:nr)
        end if
        pp = ishft(pp, -1)
        if (pp > 0_int64) then
            call z_mul_big(a, na, a, na, t, nt)
            na = nt
            a(1:na) = t(1:na)
        end if
    end do
    nout = nr
    out(1:nout) = r(1:nout)
end subroutine z_pow

integer(int64) function z_bitlen(d, n) result(L)
    integer(int64), intent(in) :: d(:)
    integer, intent(in) :: n
    integer(int64) :: v
    L = 0_int64
    if (n < 1) return
    v = d(n)
    L = int(n - 1, int64) * 16_int64
    do while (v > 0_int64)
        L = L + 1_int64
        v = ishft(v, -1)
    end do
end function z_bitlen

integer function z_bit(d, n, k) result(b)
    integer(int64), intent(in) :: d(:), k
    integer, intent(in) :: n
    integer(int64) :: limb
    integer :: off, idx
    b = 0
    if (k < 0_int64) return
    limb = k / 16_int64
    off = int(mod(k, 16_int64))
    idx = int(limb) + 1
    if (idx < 1 .or. idx > n) return
    b = int(iand(ishft(d(idx), -off), 1_int64))
end function z_bit

logical function z_sticky_below(d, n, k) result(s)
    integer(int64), intent(in) :: d(:), k
    integer, intent(in) :: n
    integer :: i, full, off, idx
    integer(int64) :: mask
    s = .false.
    if (k <= 0_int64) return
    full = int(k / 16_int64)
    do i = 1, min(full, n)
        if (d(i) /= 0_int64) then
            s = .true.
            return
        end if
    end do
    off = int(mod(k, 16_int64))
    if (off > 0) then
        idx = full + 1
        if (idx <= n) then
            mask = ishft(1_int64, off) - 1_int64
            if (iand(d(idx), mask) /= 0_int64) s = .true.
        end if
    end if
end function z_sticky_below

subroutine emit(neg, bits, y)
    logical, intent(in) :: neg
    integer(int64), intent(in) :: bits
    real(dp), intent(out) :: y
    integer(int64) :: b
    b = bits
    if (neg) b = ior(b, ishft(1_int64, 63))
    y = transfer(b, 0.0_dp)
end subroutine emit

subroutine z_round(neg, d, n, exp, y)
    logical, intent(in) :: neg
    integer(int64), intent(in) :: d(:), exp
    integer, intent(in) :: n
    real(dp), intent(out) :: y
    integer(int64) :: L, e, shift, quot, sh, drop, j
    integer :: rb
    logical :: sticky
    if (n <= 0) then
        call emit(neg, 0_int64, y)
        return
    end if
    L = z_bitlen(d, n)
    e = exp + L - 1_int64
    if (e > 1023_int64) then
        call emit(neg, ishft(2047_int64, 52), y)
        return
    end if
    if (e >= -1022_int64) then
        shift = L - 53_int64
        quot = 0_int64
        do j = 0_int64, 52_int64
            if (z_bit(d, n, L - 1_int64 - j) == 1) quot = ior(quot, ishft(1_int64, int(52_int64 - j)))
        end do
        if (shift > 0_int64) then
            rb = z_bit(d, n, L - 54_int64)
            sticky = z_sticky_below(d, n, L - 54_int64)
            if (rb == 1 .and. (sticky .or. iand(quot, 1_int64) == 1_int64)) quot = quot + 1_int64
        end if
        if (quot == ishft(1_int64, 53)) then
            quot = ishft(quot, -1)
            e = e + 1_int64
        end if
        if (e > 1023_int64) then
            call emit(neg, ishft(2047_int64, 52), y)
            return
        end if
        call emit(neg, ior(ishft(e + 1023_int64, 52), quot - ishft(1_int64, 52)), y)
        return
    end if
    sh = exp + 1074_int64
    if (sh >= 0_int64) then
        quot = 0_int64
        do j = 0_int64, L - 1_int64
            if (z_bit(d, n, j) == 1) quot = ior(quot, ishft(1_int64, int(j + sh)))
        end do
        if (quot == 0_int64) then
            call emit(neg, 0_int64, y)
        else if (quot >= ishft(1_int64, 52)) then
            call round_int(neg, quot, -1074_int64, y)
        else
            call emit(neg, quot, y)
        end if
        return
    end if
    drop = -sh
    if (drop > L + 2_int64) then
        call emit(neg, 0_int64, y)
        return
    end if
    quot = 0_int64
    do j = drop, L - 1_int64
        if (z_bit(d, n, j) == 1) quot = ior(quot, ishft(1_int64, int(j - drop)))
    end do
    rb = 0
    sticky = .false.
    if (drop >= 1_int64) rb = z_bit(d, n, drop - 1_int64)
    if (drop > 1_int64) sticky = z_sticky_below(d, n, drop - 1_int64)
    if (rb == 1 .and. (sticky .or. iand(quot, 1_int64) == 1_int64)) quot = quot + 1_int64
    if (quot == 0_int64) then
        call emit(neg, 0_int64, y)
    else if (quot >= ishft(1_int64, 52)) then
        call round_int(neg, quot, -1074_int64, y)
    else
        call emit(neg, quot, y)
    end if
end subroutine z_round

subroutine round_int(neg, numer, exp, y)
    logical, intent(in) :: neg
    integer(int64), intent(in) :: numer, exp
    real(dp), intent(out) :: y
    integer(int64) :: d(4)
    integer :: n
    call z_set_int(d, n, numer)
    call z_round(neg, d, n, exp, y)
end subroutine round_int

subroutine fp_zero(d)
    integer(int64), intent(out) :: d(:)
    d = 0_int64
end subroutine fp_zero

subroutine fp_set_one(d)
    integer(int64), intent(out) :: d(:)
    call fp_zero(d)
    d(fb / 16 + 1) = 1_int64
end subroutine fp_set_one

logical function fp_zero_p(d) result(z)
    integer(int64), intent(in) :: d(:)
    z = all(d(1:fl) == 0_int64)
end function fp_zero_p

logical function fp_low_nonzero(prod, nlow) result(z)
    integer(int64), intent(in) :: prod(:)
    integer, intent(in) :: nlow
    integer :: i
    z = .false.
    do i = 1, nlow
        if (prod(i) /= 0_int64) then
            z = .true.
            return
        end if
    end do
end function fp_low_nonzero

subroutine fp_add_small(d, m)
    integer(int64), intent(inout) :: d(:)
    integer(int64), intent(in) :: m
    integer(int64) :: acc, cur
    integer :: i
    acc = m
    i = 1
    do while (acc > 0_int64 .and. i <= fl)
        cur = d(i) + acc
        d(i) = iand(cur, zmask)
        acc = ishft(cur, -16)
        i = i + 1
    end do
end subroutine fp_add_small

subroutine fp_add(a, b, c)
    integer(int64), intent(in) :: a(:), b(:)
    integer(int64), intent(out) :: c(:)
    integer(int64) :: acc, cur
    integer :: i
    acc = 0_int64
    do i = 1, fl
        cur = a(i) + b(i) + acc
        c(i) = iand(cur, zmask)
        acc = ishft(cur, -16)
    end do
    if (acc /= 0_int64) error stop "fp_add overflow"
end subroutine fp_add

subroutine fp_shl1(d)
    integer(int64), intent(inout) :: d(:)
    integer(int64) :: carry, cur, bit
    integer :: i
    carry = 0_int64
    do i = 1, fl
        bit = ishft(d(i), -15)
        cur = ishft(d(i), 1) + carry
        d(i) = iand(cur, zmask)
        carry = bit
    end do
end subroutine fp_shl1

subroutine fp_mul_raw(a, b, prod)
    integer(int64), intent(in) :: a(:), b(:)
    integer(int64), intent(out) :: prod(:)
    integer(int64) :: acc, piece
    integer :: i, j, k
    prod = 0_int64
    do i = 1, fl
        if (a(i) == 0_int64) cycle
        acc = 0_int64
        do j = 1, fl
            piece = prod(i + j - 1) + a(i) * b(j) + acc
            prod(i + j - 1) = iand(piece, zmask)
            acc = ishft(piece, -16)
        end do
        k = i + fl
        do while (acc > 0_int64)
            if (k > size(prod)) error stop "fp_mul_raw overflow"
            piece = prod(k) + acc
            prod(k) = iand(piece, zmask)
            acc = ishft(piece, -16)
            k = k + 1
        end do
    end do
end subroutine fp_mul_raw

subroutine fp_take(prod, c, ceil)
    integer(int64), intent(in) :: prod(:)
    integer(int64), intent(out) :: c(:)
    logical, intent(in) :: ceil
    integer :: i, sh
    sh = fb / 16
    call fp_zero(c)
    do i = 1, fl
        c(i) = prod(i + sh)
    end do
    if (ceil .and. fp_low_nonzero(prod, sh)) call fp_add_small(c, 1_int64)
end subroutine fp_take

subroutine fp_div_small(a, m, c, rem)
    integer(int64), intent(in) :: a(:)
    integer(int64), intent(in) :: m
    integer(int64), intent(out) :: c(:)
    integer(int64), intent(out) :: rem
    integer(int64) :: cur
    integer :: i
    rem = 0_int64
    call fp_zero(c)
    do i = fl, 1, -1
        cur = ishft(rem, 16) + a(i)
        c(i) = cur / m
        rem = cur - c(i) * m
    end do
end subroutine fp_div_small

subroutine fp_ln2(acc)
    integer(int64), intent(out) :: acc(:)
    integer(int64) :: one(fl), s(fl), s2(fl), term(fl), nxt(fl), sum(fl), pow(fl)
    integer(int64) :: prod(2 * fl + 4), rem
    integer :: k
    call fp_set_one(one)
    call fp_div_small(one, 3_int64, s, rem)
    call fp_mul_raw(s, s, prod)
    call fp_take(prod, s2, .false.)
    pow = s
    sum = s
    do k = 1, 900
        call fp_mul_raw(pow, s2, prod)
        call fp_take(prod, pow, .false.)
        call fp_div_small(pow, int(2 * k + 1, int64), term, rem)
        if (fp_zero_p(term)) exit
        call fp_add(sum, term, nxt)
        sum = nxt
    end do
    if (.not. fp_zero_p(term)) error stop "ln2 series short"
    call fp_shl1(sum)
    acc = sum
end subroutine fp_ln2

subroutine fp_exp(z, acc)
    integer(int64), intent(in) :: z(:)
    integer(int64), intent(out) :: acc(:)
    integer(int64) :: term(fl), nxt(fl), sum(fl), prod(2 * fl + 4), rem
    integer :: k
    call fp_set_one(sum)
    term = sum
    do k = 1, 1200
        call fp_mul_raw(term, z, prod)
        call fp_take(prod, nxt, .false.)
        call fp_div_small(nxt, int(k, int64), term, rem)
        if (fp_zero_p(term)) exit
        call fp_add(sum, term, nxt)
        sum = nxt
    end do
    if (.not. fp_zero_p(term)) error stop "exp series short"
    acc = sum
end subroutine fp_exp

subroutine fp_as_z(fp, z, n)
    integer(int64), intent(in) :: fp(:)
    integer(int64), intent(out) :: z(:)
    integer, intent(out) :: n
    integer :: i
    n = fl
    do i = 1, fl
        z(i) = fp(i)
    end do
    call z_trim(z, n)
end subroutine fp_as_z

logical function fp_ge(a, b) result(ge)
    integer(int64), intent(in) :: a(:), b(:)
    integer :: i
    ge = .true.
    do i = fl, 1, -1
        if (a(i) > b(i)) return
        if (a(i) < b(i)) then
            ge = .false.
            return
        end if
    end do
end function fp_ge

subroutine fp_sub(a, b)
    integer(int64), intent(inout) :: a(:)
    integer(int64), intent(in) :: b(:)
    integer(int64) :: cur
    integer :: i, borrow
    borrow = 0
    do i = 1, fl
        cur = a(i) - b(i) - borrow
        if (cur < 0_int64) then
            cur = cur + 65536_int64
            borrow = 1
        else
            borrow = 0
        end if
        a(i) = cur
    end do
end subroutine fp_sub

subroutine fp_recip_unit(den, r)
    ! r = floor(2**(2*fb) / den). den is in about [2**fb, 2**(fb+1)).
    integer(int64), intent(in) :: den(:)
    integer(int64), intent(out) :: r(:)
    integer(int64) :: rem(fl)
    integer :: bit, limb, off
    call fp_zero(rem)
    call fp_zero(r)
    do bit = 2 * fb, 0, -1
        call fp_shl1(rem)
        if (bit == 2 * fb) call fp_add_small(rem, 1_int64)
        if (fp_ge(rem, den)) then
            call fp_sub(rem, den)
            limb = bit / 16 + 1
            off = mod(bit, 16)
            if (limb <= fl) r(limb) = ior(r(limb), ishft(1_int64, off))
        end if
    end do
end subroutine fp_recip_unit

subroutine fp_div_i64(a, den, q)
    integer(int64), intent(in) :: a(:), den
    integer(int64), intent(out) :: q(:)
    integer(int64) :: rem, cur, qh, ql
    integer :: i
    rem = 0_int64
    call fp_zero(q)
    do i = fl, 1, -1
        cur = ishft(rem, 8) + ishft(a(i), -8)
        qh = cur / den
        rem = cur - qh * den
        cur = ishft(rem, 8) + iand(a(i), 255_int64)
        ql = cur / den
        rem = cur - ql * den
        q(i) = ior(ishft(qh, 8), ql)
    end do
end subroutine fp_div_i64

subroutine fp_mul_i64(a, m, c)
    integer(int64), intent(in) :: a(:), m
    integer(int64), intent(out) :: c(:)
    integer(int64) :: acc, piece, parts(4)
    integer :: s, i
    call fp_zero(c)
    parts(1) = iand(m, zmask)
    parts(2) = iand(ishft(m, -16), zmask)
    parts(3) = iand(ishft(m, -32), zmask)
    parts(4) = iand(ishft(m, -48), zmask)
    do s = 1, 4
        if (parts(s) == 0_int64) cycle
        acc = 0_int64
        do i = 1, fl
            piece = a(i) * parts(s) + acc
            if (i + s - 1 > fl) then
                if (piece /= 0_int64) error stop "fp_mul_i64 overflow"
                exit
            end if
            piece = piece + c(i + s - 1)
            c(i + s - 1) = iand(piece, zmask)
            acc = ishft(piece, -16)
        end do
        if (acc /= 0_int64) error stop "fp_mul_i64 carry"
    end do
end subroutine fp_mul_i64

subroutine fp_div_i64_small(a, m, c)
    integer(int64), intent(in) :: a(:)
    integer(int64), intent(in) :: m
    integer(int64), intent(out) :: c(:)
    integer(int64) :: rem
    call fp_div_small(a, m, c, rem)
end subroutine fp_div_i64_small

integer function z_bitlen_small(v) result(L)
    integer(int64), intent(in) :: v
    integer(int64) :: x
    L = 0
    x = v
    do while (x > 0_int64)
        L = L + 1
        x = ishft(x, -1)
    end do
end function z_bitlen_small

subroutine place_i64_at(d, v, sh)
    integer(int64), intent(inout) :: d(:)
    integer(int64), intent(in) :: v
    integer, intent(in) :: sh
    integer :: i, limb, off
    integer(int64) :: x
    x = v
    i = 0
    do while (x > 0_int64)
        limb = (sh + i) / 16 + 1
        off = mod(sh + i, 16)
        if (limb > fl) error stop "place overflow"
        d(limb) = ior(d(limb), ishft(iand(x, 1_int64), off))
        ! Low bit of v is placed at bit sh.
        x = ishft(x, -1)
        i = i + 1
    end do
end subroutine place_i64_at

subroutine split_scaled(qv, rshift, ia, frac, pure_int, overflow)
    integer(int64), intent(in) :: qv(:)
    integer, intent(in) :: rshift
    integer(int64), intent(out) :: ia
    integer(int64), intent(out) :: frac(:)
    logical, intent(out) :: pure_int, overflow
    integer :: i, bit, topb
    pure_int = .false.
    overflow = .false.
    call fp_zero(frac)
    ia = 0_int64
    topb = -1
    do i = fl, 1, -1
        if (qv(i) /= 0_int64) then
            topb = (i - 1) * 16 + z_bitlen_small(qv(i)) - 1
            exit
        end if
    end do
    if (topb < 0) then
        pure_int = .true.
        ia = 0_int64
        return
    end if
    ! real = qv / 2**(-rshift) if rshift negative?
    ! |y| = qv_int * 2**rshift, rshift = pe-fb usually negative.
    if (rshift >= 0) then
        if (topb + rshift > 40) then
            overflow = .true.
            return
        end if
        pure_int = .true.
        ia = 0_int64
        ! small left shift, read integer
        ia = read_shifted(qv, -rshift)
        return
    end if
    ! right shift by -rshift
    bit = -rshift
    if (topb < bit) then
        ! |y| < 1
        ia = 0_int64
        call frac_from_shift(qv, bit, frac)
        if (fp_zero_p(frac)) pure_int = .true.
        return
    end if
    if (topb - bit > 40) then
        overflow = .true.
        return
    end if
    ia = read_shifted(qv, bit)
    call frac_from_shift(qv, bit, frac)
    if (fp_zero_p(frac)) pure_int = .true.
end subroutine split_scaled

integer(int64) function read_shifted(d, sh) result(v)
    ! floor(int(d) / 2**sh) if sh>=0, or int(d)<<(-sh) if sh<0. Result must fit int64.
    integer(int64), intent(in) :: d(:)
    integer, intent(in) :: sh
    integer :: bit, limb, off, nbits
    v = 0_int64
    if (sh >= 0) then
        nbits = 0
        do bit = sh, sh + 62
            limb = bit / 16 + 1
            off = mod(bit, 16)
            if (limb <= fl) then
                if (iand(ishft(d(limb), -off), 1_int64) == 1_int64) v = ior(v, ishft(1_int64, nbits))
            end if
            nbits = nbits + 1
        end do
    else
        error stop "read_shifted left"
    end if
end function read_shifted

subroutine frac_from_shift(d, sh, frac)
    ! fractional part of int(d)/2**sh as fixed point at binary point fb.
    integer(int64), intent(in) :: d(:)
    integer, intent(in) :: sh
    integer(int64), intent(out) :: frac(:)
    integer :: bit, limb, off, dest
    call fp_zero(frac)
    if (sh <= 0) return
    do bit = 0, fb - 1
        dest = sh - 1 - bit
        if (dest < 0) exit
        limb = dest / 16 + 1
        off = mod(dest, 16)
        if (limb <= fl) then
            if (iand(ishft(d(limb), -off), 1_int64) == 1_int64) then
                frac((fb - 1 - bit) / 16 + 1) = ior(frac((fb - 1 - bit) / 16 + 1), ishft(1_int64, mod(fb - 1 - bit, 16)))
            end if
        end if
    end do
end subroutine frac_from_shift

subroutine scale_clamp_mant(neg, wm, we, nfloor, y)
    logical, intent(in) :: neg
    integer(int64), intent(in) :: wm, we, nfloor
    real(dp), intent(out) :: y
    integer(int64) :: d(8)
    integer :: n
    call z_set_int(d, n, wm)
    call z_round(neg, d, n, we + nfloor, y)
end subroutine scale_clamp_mant
subroutine factor_hpow(p, half, status, nfloor, sig)
    real(dp), intent(in) :: p, half
    integer, intent(out) :: status
    integer(int64), intent(out) :: nfloor
    integer(int64), intent(out) :: sig(:)
    integer(int64) :: hm, he, num, den, pm, pe
    integer(int64) :: s(fl), s2(fl), pow(fl), term(fl), nxt(fl), lsig(fl)
    integer(int64) :: l2(fl), inv(fl), lg(fl), mag(fl), prod(2 * fl + 4)
    integer(int64) :: qv(fl), frac(fl), ex(fl)
    integer(int64) :: eun, ia
    integer :: k
    logical :: negy, negsum, pure_int, overflow
    status = 0
    nfloor = 0_int64
    call fp_zero(sig)
    if (half == 0.0_dp) then
        ! Flushed half-width is the positive gap below the tiniest subnormal.
        hm = ishft(1_int64, 52)
        he = -1075_int64 - 52_int64
    else
        call split_dp(half, hm, he)
    end if
    do while (hm > 0_int64 .and. hm < ishft(1_int64, 52))
        hm = ishft(hm, 1)
        he = he - 1_int64
    end do
    eun = he + 52_int64
    num = hm - ishft(1_int64, 52)
    den = hm + ishft(1_int64, 52)
    ! numerator num << fb
    if (num == 0_int64) then
        call fp_zero(s)
    else
        call fp_zero(qv)
        qv(1) = iand(num, zmask)
        qv(2) = iand(ishft(num, -16), zmask)
        qv(3) = iand(ishft(num, -32), zmask)
        qv(4) = iand(ishft(num, -48), zmask)
        do k = 1, fb
            call fp_shl1(qv)
        end do
        call fp_div_i64(qv, den, s)
    end if
    call fp_mul_raw(s, s, prod)
    call fp_take(prod, s2, .false.)
    pow = s
    lsig = s
    do k = 1, 900
        call fp_mul_raw(pow, s2, prod)
        call fp_take(prod, pow, .false.)
        call fp_div_i64_small(pow, int(2 * k + 1, int64), term)
        if (fp_zero_p(term)) exit
        call fp_add(lsig, term, nxt)
        lsig = nxt
    end do
    call fp_shl1(lsig)
    call fp_ln2(l2)
    call fp_recip_unit(l2, inv)
    call fp_mul_raw(lsig, inv, prod)
    call fp_take(prod, lg, .false.)
    ! |sum| = |eun + log2(sig)| with log2 in lg
    if (eun >= 0_int64) then
        negsum = .false.
        call fp_zero(mag)
        if (eun /= 0_int64) then
            call place_i64_at(mag, eun, fb)
        end if
        call fp_add(mag, lg, nxt)
        mag = nxt
    else
        negsum = .true.
        call fp_zero(mag)
        call place_i64_at(mag, -eun, fb)
        call fp_sub(mag, lg)
        if (fp_zero_p(mag)) negsum = .false.
    end if
    call split_dp(abs(p), pm, pe)
    negy = negsum .neqv. (p < 0.0_dp)
    call fp_mul_i64(mag, pm, qv)
    ! |y| = qv * 2**(pe - fb)
    pure_int = .false.
    call split_scaled(qv, int(pe - int(fb, int64)), ia, frac, pure_int, overflow)
    if (overflow) then
        if (negy) then
            status = 3
        else
            status = 2
        end if
        return
    end if
    if (pure_int) then
        if (negy) then
            nfloor = -ia
        else
            nfloor = ia
        end if
        status = 1
        return
    end if
    call fp_ln2(l2)
    call fp_mul_raw(frac, l2, prod)
    call fp_take(prod, term, .false.)
    call fp_exp(term, ex)
    if (negy) then
        call fp_recip_unit(ex, inv)
        nfloor = -ia
        sig = inv
    else
        nfloor = ia
        sig = ex
    end if
    status = 0
end subroutine factor_hpow

subroutine apply_sig(w, status, nfloor, sig, y)
    real(dp), intent(in) :: w
    integer, intent(in) :: status
    integer(int64), intent(in) :: nfloor
    integer(int64), intent(in) :: sig(:)
    real(dp), intent(out) :: y
    integer(int64) :: wm, we, dust(fl), plo(zmax), phi(zmax)
    integer :: nplo, nphi
    logical :: negw
    real(dp) :: ylo, yhi
    negw = w < 0.0_dp
    if (w == 0.0_dp) then
        y = w
        return
    end if
    if (status == 2) then
        call emit(negw, ishft(2047_int64, 52), y)
        return
    end if
    if (status == 3) then
        call emit(negw, 0_int64, y)
        return
    end if
    call split_dp(abs(w), wm, we)
    if (status == 1) then
        call scale_clamp_mant(negw, wm, we, nfloor, y)
        return
    end if
    dust = sig
    dust(3) = ior(dust(3), ishft(1_int64, 8))
    call fp_as_z(sig, plo, nplo)
    call z_mul_int(plo, nplo, wm)
    call fp_as_z(dust, phi, nphi)
    call z_mul_int(phi, nphi, wm)
    call z_round(negw, plo, nplo, we + nfloor - int(fb, int64), ylo)
    call z_round(negw, phi, nphi, we + nfloor - int(fb, int64), yhi)
    if (transfer(ylo, 0_int64) /= transfer(yhi, 0_int64)) error stop "apply amb"
    y = ylo
end subroutine apply_sig

!> True when odd(H)**p fits in the limb buffer.
logical function int_pow_fits(half, p) result(ok)
    real(dp), intent(in) :: half, p
    integer(int64) :: hm, he, bits, x
    if (half == 0.0_dp) then
        hm = 1_int64
    else
        call split_dp(half, hm, he)
    end if
    do while (hm > 1_int64 .and. iand(hm, 1_int64) == 0_int64)
        hm = ishft(hm, -1)
    end do
    bits = 0_int64
    x = hm
    do while (x > 0_int64)
        bits = bits + 1_int64
        x = ishft(x, -1)
    end do
    ok = bits * int(p, int64) + 256_int64 < int(zmax, int64) * 16_int64
end function int_pow_fits

!> Exact odd significand of H**p and the power of two that goes with it.
subroutine factor_int(p, half, fexp, mag, nmag)
    integer(int64), intent(in) :: p
    real(dp), intent(in) :: half
    integer(int64), intent(out) :: fexp
    integer(int64), intent(out) :: mag(:)
    integer, intent(out) :: nmag
    integer(int64) :: hm, he
    if (half == 0.0_dp) then
        hm = 1_int64
        he = -1075_int64
    else
        call split_dp(half, hm, he)
    end if
    do while (hm > 0_int64 .and. iand(hm, 1_int64) == 0_int64)
        hm = ishft(hm, -1)
        he = he + 1_int64
    end do
    call z_pow(hm, p, mag, nmag)
    fexp = he * p
end subroutine factor_int

!> One ties-to-even rounding of w times an exact integer power.
subroutine apply_int(w, fexp, mag, nmag, y)
    real(dp), intent(in) :: w
    integer(int64), intent(in) :: fexp
    integer(int64), intent(in) :: mag(:)
    integer, intent(in) :: nmag
    real(dp), intent(out) :: y
    integer(int64) :: d(zmax), wm, we
    integer :: nd
    logical :: neg
    if (w == 0.0_dp) then
        y = w
        return
    end if
    neg = w < 0.0_dp
    call split_dp(abs(w), wm, we)
    nd = nmag
    if (nd > 0) d(1:nd) = mag(1:nd)
    call z_mul_int(d, nd, wm)
    call z_round(neg, d, nd, fexp + we, y)
end subroutine apply_int

!> Correctly rounded w*H**p. The heavy factor is prepared once by the caller.
!> mode 1 is scale() with the historical ±4096 clamp (power-of-two H, integer p).
!> mode 2 is one IEEE multiply (p = 1). mode 4 is an exact integer power.
!> mode 5 rounds a fixed-point factor. No extended real kind.
!> -ffast-math and -ffinite-math-only are unsupported: they delete
!> ieee_is_finite, so a non-finite node or weight is not rejected.
subroutine scaled_weight(w, p, half, mode, iscal, fac, sigst, signf, wsig, fexp, mag, nmag, tw)
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
    real(dp), intent(in) :: w, p, half, fac
    integer, intent(in) :: mode, iscal, sigst, nmag
    integer(int64), intent(in) :: signf, fexp
    integer(int64), intent(in) :: wsig(:)
    integer(int64), intent(in) :: mag(:)
    real(dp), intent(out) :: tw
    if (w == 0.0_dp .or. p == 0.0_dp .or. half == 1.0_dp) then
        tw = w
        return
    end if
    select case (mode)
    case (1)
        tw = scale(w, iscal)
    case (2)
        tw = fac * w
    case (3)
        if (half > 1.0_dp) then
            tw = ieee_value(0.0_dp, ieee_positive_inf)
        else
            tw = 0.0_dp
        end if
        if (w < 0.0_dp) tw = -tw
    case (4)
        call apply_int(w, fexp, mag, nmag, tw)
    case (5)
        call apply_sig(w, sigst, signf, wsig, tw)
    case default
        tw = w
    end select
end subroutine scaled_weight

end module GaussJacobiQuad

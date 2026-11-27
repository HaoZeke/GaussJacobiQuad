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
implicit none
public :: gauss_jacobi, gauss_jacobi_rule, gauss_jacobi_batch_caf
public :: select_method_auto, gauss_jacobi_scale
private :: scale_interval, gjp_scale_f, affine_image, weight_image
private :: aff_from_double, aff_uadd, aff_usub, aff_ucmp, aff_sadd, aff_smul
private :: aff_bit, aff_sticky, aff_round, pow2_exp
! Each operand is an integer significand times 2**1074. aff_round divides
! that product by 2**2149. 280 base-2**16 limbs cover every finite binary64
! product (the leading limb stays near 262).
integer, parameter, private :: aff_nl = 280
integer, parameter, private :: aff_lw = 16
integer(int64), parameter, private :: aff_mask = 65535_int64
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
    real(dp) :: half, mid, p, t, tw
    real(dp), allocatable :: xt(:), wt(:)
    integer :: i
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
        tw = weight_image(wts(i), p, half)
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

!> Correctly rounded w * H**p. H is half, or 2**-1075 when half flushed to 0.
!> A power-of-two H and an integer p use scale, so overflow is exact.
!> Every other product is formed in real128 and rounded to binary64.
!> A non-finite result is an overflow; the caller does not write.
pure function weight_image(w, p, half) result(tw)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_positive_inf
    use, intrinsic :: iso_fortran_env, only: real128
    real(dp), intent(in) :: w, p, half
    real(dp) :: tw
    real(real128) :: hq, wq, pq, prod
    integer :: k, iscal
    integer(int64) :: pk
    ! Factor 1 skips the product, including a non-finite power times one.
    if (w == 0.0_dp .or. p == 0.0_dp .or. half == 1.0_dp) then
        tw = w
        return
    end if
    if (.not. ieee_is_finite(p)) then
        if (half > 1.0_dp) then
            tw = ieee_value(0.0_dp, ieee_positive_inf)
        else
            tw = 0.0_dp
        end if
        if (w < 0.0_dp) tw = -tw
        return
    end if
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
        tw = scale(w, iscal)
        return
    end if
    if (half == 0.0_dp) then
        hq = 2.0_real128**(-1075)
    else
        hq = real(half, real128)
    end if
    wq = abs(real(w, real128))
    pq = real(p, real128)
    prod = exp(log(wq) + pq * log(hq))
    tw = real(prod, dp)
    if (w < 0.0_dp) tw = -tw
end function weight_image

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

end module GaussJacobiQuad

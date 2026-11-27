! BEGIN_HEADER
! -----------------------------------------------------------------------------
! Gauss-Jacobi Quadrature Implementation
! Authors: Rohit Goswami <rgoswami[at]ieee.org>
! Source: GaussJacobiQuad Library
! License: MIT
! GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
! -----------------------------------------------------------------------------
! END_HEADER
!> Integrate a polynomial on a non-symmetric interval with the scaled rule.
program test_scale_interval
use GaussJacobiQuad, only: gauss_jacobi_rule, gauss_jacobi_scale
use gjp_types, only: dp
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_next_after
use, intrinsic :: iso_c_binding, only: c_double, c_int
use, intrinsic :: iso_fortran_env, only: int64
implicit none

interface
    integer(c_int) function gjp_scale_f(npts, alpha, beta, a, b, x, wts) &
        bind(C, name="gjp_scale_f")
        import :: c_int, c_double
        implicit none
        integer(c_int), value :: npts
        real(c_double), value :: alpha, beta, a, b
        real(c_double) :: x(npts), wts(npts)
    end function gjp_scale_f
end interface

integer, parameter :: n = 8
integer, parameter :: n4 = 4
real(dp), parameter :: alpha = 1.0_dp, beta = 2.0_dp
real(dp), parameter :: a = 0.0_dp, b = 3.0_dp
real(dp) :: x(n), w(n), quad, exact, err, wsum
real(dp) :: x4(n4), w4(n4), y4(n4), z4(n4), lo, hi
real(dp) :: xs(1), ws(1), xkeep, wkeep
integer :: i
integer(c_int) :: stc
real(dp) :: xu(4), wu(4), xl(8), wl(8), edge

call gauss_jacobi_rule(n, alpha, beta, x, w, "rec")
call gauss_jacobi_scale(n, alpha, beta, a, b, x, w)

quad = 0.0_dp
do i = 1, n
    quad = quad + w(i) * x(i)**3
end do
! ∫_0^3 (3-t)^1 t^2 t^3 dt = ∫_0^3 (3-t) t^5 dt = 729/14
exact = 729.0_dp / 14.0_dp
err = abs(quad - exact)
! f = 1 moment: (b-a)^(alpha+beta+1) * B(beta+1, alpha+1) = 81 * B(3,2) = 27/4
wsum = sum(w)
print '(A,ES16.8,A,ES16.8)', "scale [0,3] poly err=", err, " moment err=", abs(wsum - 27.0_dp / 4.0_dp)
if (err > 1.0e-12_dp) error stop "scaled polynomial integral failed"
if (abs(wsum - 27.0_dp / 4.0_dp) > 1.0e-12_dp) error stop "scaled weight sum failed"
if (any(x <= a) .or. any(x >= b)) error stop "scaled nodes outside (a,b)"
if (any(w <= 0.0_dp)) error stop "scaled weights not positive"

! Huge finite endpoints: finite nodes and weights inside (a, b).
! Status is not returned here; a rejected scale error-stops.
call gauss_jacobi_rule(n4, 0.0_dp, 0.0_dp, x4, w4, "rec")
y4 = x4
z4 = w4
lo = 1.0e308_dp
hi = 1.7e308_dp
call gauss_jacobi_scale(n4, 0.0_dp, 0.0_dp, lo, hi, x4, w4)
if (.not. all(ieee_is_finite(x4)) .or. .not. all(ieee_is_finite(w4))) &
    error stop "non-finite values on [1e308, 1.7e308]"
if (any(x4 <= lo) .or. any(x4 >= hi)) error stop "nodes outside [1e308, 1.7e308]"
lo = -1.0e308_dp
hi = 1.0e308_dp
call gauss_jacobi_scale(n4, 0.0_dp, 0.0_dp, lo, hi, y4, z4)
if (.not. all(ieee_is_finite(y4)) .or. .not. all(ieee_is_finite(z4))) &
    error stop "non-finite values on [-1e308, 1e308]"
if (any(y4 <= lo) .or. any(y4 >= hi)) error stop "nodes outside [-1e308, 1e308]"

! Endpoints that the lease sum sends outside [a, b] or to infinity.
xs(1) = -1.0_dp
ws(1) = 1.0_dp
stc = gjp_scale_f(1_c_int, 0.0_c_double, 0.0_c_double, 1.0_c_double, 1.0e16_c_double, xs, ws)
if (stc /= 0_c_int .or. .not. ieee_is_finite(xs(1)) .or. xs(1) /= 1.0_dp) &
    error stop "node [1, 1e16] x=-1"
print '(A,ES24.16)', "endpoint [1,1e16] x=-1 ", xs(1)

xs(1) = ieee_next_after(-1.0_dp, 0.0_dp)
ws(1) = 1.0_dp
stc = gjp_scale_f(1_c_int, 0.0_c_double, 0.0_c_double, 1.0_c_double, 1.0e16_c_double, xs, ws)
if (stc /= 0_c_int .or. .not. ieee_is_finite(xs(1))) error stop "node [1, 1e16] nextafter status"
if (xs(1) <= 1.0_dp .or. xs(1) >= 1.0e16_dp) error stop "node [1, 1e16] nextafter outside"
print '(A,ES24.16)', "endpoint [1,1e16] next ", xs(1)

xs(1) = -1.0_dp
ws(1) = 1.0_dp
stc = gjp_scale_f(1_c_int, 0.0_c_double, 0.0_c_double, 1.0e6_c_double, 1.0e20_c_double, xs, ws)
if (stc /= 0_c_int .or. .not. ieee_is_finite(xs(1)) .or. xs(1) /= 1.0e6_dp) &
    error stop "node [1e6, 1e20] x=-1"
print '(A,ES24.16)', "endpoint [1e6,1e20] x=-1 ", xs(1)

xs(1) = 1.0_dp
ws(1) = 1.0_dp
stc = gjp_scale_f(1_c_int, 0.0_c_double, 0.0_c_double, 1.0e308_c_double, huge(1.0_dp), xs, ws)
if (stc /= 0_c_int .or. .not. ieee_is_finite(xs(1)) .or. xs(1) /= huge(1.0_dp)) &
    error stop "node [1e308, huge] x=1"
print '(A,ES24.16)', "endpoint [1e308,huge] x=1 ", xs(1)

! Finite factor 1e308 times weight 2 overflows. Status 7, buffers unchanged.
xs(1) = 0.0_dp
ws(1) = 2.0_dp
xkeep = xs(1)
wkeep = ws(1)
stc = gjp_scale_f(1_c_int, 0.0_c_double, 0.0_c_double, -1.0e308_c_double, 1.0e308_c_double, xs, ws)
if (stc /= 7_c_int) error stop "weight overflow status"
if (xs(1) /= xkeep .or. ws(1) /= wkeep) error stop "weight overflow wrote buffers"
print '(A,I0)', "weight overflow status ", stc

! Endpoints the non-fused lease leaves strictly inside (a, b).
xs(1) = -1.0_dp
ws(1) = 0.0_dp
call gauss_jacobi_scale(1, 0.0_dp, 0.0_dp, -1.0_dp, 1.0e16_dp, xs, ws)
if (xs(1) /= -1.0_dp) error stop "node [-1, 1e16] x=-1"
xs(1) = ieee_next_after(-1.0_dp, 0.0_dp)
ws(1) = 0.0_dp
call gauss_jacobi_scale(1, 0.0_dp, 0.0_dp, -1.0_dp, 1.0e16_dp, xs, ws)
if (xs(1) >= 0.0_dp) error stop "node [-1, 1e16] nextafter sign"
if (abs(xs(1) + 0.44488848768742173_dp) > 1.0e-15_dp) error stop "node [-1, 1e16] nextafter"

xs(1) = 1.0_dp
ws(1) = 0.0_dp
call gauss_jacobi_scale(1, 0.0_dp, 0.0_dp, -1.0e16_dp, 1.0_dp, xs, ws)
if (xs(1) /= 1.0_dp) error stop "node [-1e16, 1] x=1"
xs(1) = ieee_next_after(1.0_dp, -1.0_dp)
ws(1) = 0.0_dp
call gauss_jacobi_scale(1, 0.0_dp, 0.0_dp, -1.0e16_dp, 1.0_dp, xs, ws)
if (xs(1) <= 0.0_dp) error stop "node [-1e16, 1] nextafter sign"
if (abs(xs(1) - 0.44488848768742173_dp) > 1.0e-15_dp) error stop "node [-1e16, 1] nextafter"

xs(1) = -1.0_dp
ws(1) = 0.0_dp
call gauss_jacobi_scale(1, 0.0_dp, 0.0_dp, -1.0e16_dp, -1.0_dp, xs, ws)
if (xs(1) /= -1.0e16_dp) error stop "node [-1e16, -1] x=-1"
xs(1) = 1.0_dp
ws(1) = 0.0_dp
call gauss_jacobi_scale(1, 0.0_dp, 0.0_dp, -1.0e16_dp, -1.0_dp, xs, ws)
if (xs(1) /= -1.0_dp) error stop "node [-1e16, -1] x=1"

xs(1) = ieee_next_after(-1.0_dp, 0.0_dp)
ws(1) = 0.0_dp
call gauss_jacobi_scale(1, 0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp, xs, ws)
if (xs(1) /= 1.6653345369377348e-16_dp) error stop "node [0, 3] nextafter"

xs(1) = 1.0_dp
ws(1) = 0.0_dp
call gauss_jacobi_scale(1, 0.0_dp, 0.0_dp, -1.0e308_dp, 1.0_dp, xs, ws)
if (xs(1) /= 1.0_dp) error stop "node [-1e308, 1] x=1"

edge = 2.0_dp**17
xs(1) = 0.3_dp
ws(1) = 0.0_dp
call gauss_jacobi_scale(1, 0.0_dp, 0.0_dp, edge, ieee_next_after(edge, huge(edge)), xs, ws)
if (xs(1) /= edge .and. xs(1) /= ieee_next_after(edge, huge(edge))) &
    error stop "node one-ulp interval"

! half**(alpha+beta+1) must not be range-checked on its own.
call gauss_jacobi_rule(4, 1022.0_dp, 0.0_dp, xu, wu, "rec")
call gauss_jacobi_scale(4, 1022.0_dp, 0.0_dp, 0.0_dp, 0.8_dp, xu, wu)
if (.not. all(ieee_is_finite(wu)) .or. any(wu == 0.0_dp)) error stop "flushed small weight"
if (any(wu < 1.0e-120_dp) .or. any(wu > 1.0e-80_dp)) error stop "small weight magnitude"

call gauss_jacobi_rule(8, 0.0_dp, 0.0_dp, xl, wl, "rec")
call gauss_jacobi_scale(8, 1023.0_dp, 0.0_dp, -2.0_dp, 2.0_dp, xl, wl)
if (.not. all(ieee_is_finite(wl)) .or. any(wl == 0.0_dp)) error stop "false weight overflow"
if (any(wl < 1.0e300_dp)) error stop "power-of-two weight magnitude"

xs(1) = 0.0_dp
ws(1) = 1.0_dp
edge = ieee_next_after(0.0_dp, 1.0_dp)
call gauss_jacobi_scale(1, -0.25_dp, -0.25_dp, 0.0_dp, edge, xs, ws)
if (.not. ieee_is_finite(ws(1)) .or. ws(1) == 0.0_dp) error stop "subnormal half p=0.5"
if (ws(1) < 1.0e-170_dp .or. ws(1) > 1.0e-150_dp) error stop "subnormal half p=0.5 magnitude"
xs(1) = 0.0_dp
ws(1) = 1.0_dp
call gauss_jacobi_scale(1, -0.9_dp, -0.9_dp, 0.0_dp, edge, xs, ws)
if (.not. ieee_is_finite(ws(1)) .or. ws(1) < 1.0e250_dp) error stop "subnormal half p=-0.8"

! Reviewer images, plus near-cancellation, overflow-adjacent, and huge-ratio
! pairs. Each entry is a, b, x, and the correctly rounded image.
block
    integer, parameter :: nlock = 14
    integer(int64), parameter :: lock(4, nlock) = reshape([ &
                                                          ior(ishft(int(z'C00B0B8A', int64), 32), int(z'D6203548', int64)), &
                                                          ior(ishft(int(z'40417FD9', int64), 32), int(z'C31752FE', int64)), &
                                                          ior(ishft(int(z'BFEFFFFF', int64), 32), int(z'FFFFFFFF', int64)), &
                                                          ior(ishft(int(z'C00B0B8A', int64), 32), int(z'D6203543', int64)), &
                                                          ior(ishft(int(z'BFA95FA5', int64), 32), int(z'F425D75E', int64)), &
                                                          ior(ishft(int(z'3FE5E952', int64), 32), int(z'758AE35F', int64)), &
                                                          ior(ishft(int(z'BFEBB093', int64), 32), int(z'A01F6698', int64)), &
                                                          ior(ishft(int(z'BF1B94A3', int64), 32), int(z'3920DC3A', int64)), &
                                                          ior(ishft(int(z'BF343DCC', int64), 32), int(z'3E402080', int64)), &
                                                          ior(ishft(int(z'3FB00880', int64), 32), int(z'ED342634', int64)), &
                                                          ior(ishft(int(z'BFEFFFFF', int64), 32), int(z'FFFFFFFF', int64)), &
                                                          ior(ishft(int(z'BF343DCC', int64), 32), int(z'3E402040', int64)), &
                                                          ior(ishft(int(z'7E37E43C', int64), 32), int(z'8800759C', int64)), &
                                                          ior(ishft(int(z'7FEFFFFF', int64), 32), int(z'FFFFFFFF', int64)), &
                                                          ior(ishft(int(z'BFEFFFFF', int64), 32), int(z'FFFFFFFF', int64)), &
                                                          ior(ishft(int(z'7E37E43C', int64), 32), int(z'8C00759C', int64)), &
                                                          ior(ishft(int(z'C1E00000', int64), 32), int(z'00000000', int64)), &
                                                          ior(ishft(int(z'45400000', int64), 32), int(z'00000000', int64)), &
                                                          ior(ishft(int(z'BFEFFFFF', int64), 32), int(z'FFFFFFFF', int64)), &
                                                          ior(ishft(int(z'3E800000', int64), 32), int(z'00000000', int64)), &
                                                          ior(ishft(int(z'00000000', int64), 32), int(z'00000000', int64)), &
                                                          ior(ishft(int(z'40080000', int64), 32), int(z'00000000', int64)), &
                                                          ior(ishft(int(z'BFEFFFFF', int64), 32), int(z'FFFFFFFF', int64)), &
                                                          ior(ishft(int(z'3CA80000', int64), 32), int(z'00000000', int64)), &
                                                          ior(ishft(int(z'BFF00000', int64), 32), int(z'00000000', int64)), &
                                                          ior(ishft(int(z'4341C379', int64), 32), int(z'37E08000', int64)), &
                                                          ior(ishft(int(z'BFEFFFFF', int64), 32), int(z'FFFFFFFF', int64)), &
                                                          ior(ishft(int(z'BFDC790D', int64), 32), int(z'903EFFFF', int64)), &
                                                          ior(ishft(int(z'4010322F', int64), 32), int(z'E7216CB2', int64)), &
                                                          ior(ishft(int(z'4010322F', int64), 32), int(z'E7216CB4', int64)), &
                                                          ior(ishft(int(z'3FD00000', int64), 32), int(z'00000000', int64)), &
                                                          ior(ishft(int(z'4010322F', int64), 32), int(z'E7216CB3', int64)), &
                                                          ior(ishft(int(z'C0F00000', int64), 32), int(z'00000000', int64)), &
                                                          ior(ishft(int(z'45200000', int64), 32), int(z'00000000', int64)), &
                                                          ior(ishft(int(z'BFEFFFFF', int64), 32), int(z'FFFFFFFF', int64)), &
                                                          ior(ishft(int(z'41BFFF00', int64), 32), int(z'00000000', int64)), &
                                                          ior(ishft(int(z'C11CBA6B', int64), 32), int(z'387D1AE4', int64)), &
                                                          ior(ishft(int(z'412A81F5', int64), 32), int(z'077331B5', int64)), &
                                                          ior(ishft(int(z'3FD160B5', int64), 32), int(z'1B47533C', int64)), &
                                                          ior(ishft(int(z'41173DC9', int64), 32), int(z'2735AF53', int64)), &
                                                          ior(ishft(int(z'C0180029', int64), 32), int(z'94F18826', int64)), &
                                                          ior(ishft(int(z'C0180029', int64), 32), int(z'94F18822', int64)), &
                                                          ior(ishft(int(z'3FD00000', int64), 32), int(z'00000000', int64)), &
                                                          ior(ishft(int(z'C0180029', int64), 32), int(z'94F18824', int64)), &
                                                          ior(ishft(int(z'C2400000', int64), 32), int(z'00000000', int64)), &
                                                          ior(ishft(int(z'43D00000', int64), 32), int(z'00000000', int64)), &
                                                          ior(ishft(int(z'BFEFFFFF', int64), 32), int(z'FFFFFFFF', int64)), &
                                                          ior(ishft(int(z'C23FFFFF', int64), 32), int(z'FF000000', int64)), &
                                                          ior(ishft(int(z'7FE1CCF3', int64), 32), int(z'85EBC8A0', int64)), &
                                                          ior(ishft(int(z'7FEFFFFF', int64), 32), int(z'FFFFFFFF', int64)), &
                                                          ior(ishft(int(z'BFEFFFFF', int64), 32), int(z'FFFFFFFF', int64)), &
                                                          ior(ishft(int(z'7FE1CCF3', int64), 32), int(z'85EBC8A0', int64)), &
                                                          ior(ishft(int(z'C12C80B9', int64), 32), int(z'7F64820E', int64)), &
                                                          ior(ishft(int(z'41237FE0', int64), 32), int(z'AC94AC8B', int64)), &
                                                          ior(ishft(int(z'BFDEFB68', int64), 32), int(z'CF907B80', int64)), &
                                                          ior(ishft(int(z'C1201ED9', int64), 32), int(z'085CF85A', int64)) &
                                                          ], [4, nlock])
    integer :: k
    do k = 1, nlock
        xs(1) = transfer(lock(3, k), 0.0_dp)
        ws(1) = 0.0_dp
        call gauss_jacobi_scale(1, 0.0_dp, 0.0_dp, &
                                transfer(lock(1, k), 0.0_dp), transfer(lock(2, k), 0.0_dp), xs, ws)
        if (transfer(xs(1), 0_int64) /= lock(4, k)) error stop "affine image"
    end do
end block
xs(1) = 0.0_dp
ws(1) = 2.0_dp
call gauss_jacobi_scale(1, huge(1.0_dp), huge(1.0_dp), 0.0_dp, 2.0_dp, xs, ws)
if (ws(1) /= 2.0_dp) error stop "unit half weight"
stc = gjp_scale_f(1_c_int, huge(0.0_c_double), huge(0.0_c_double), &
                  0.0_c_double, 2.0_c_double, xs, ws)
if (stc /= 0_c_int .or. ws(1) /= 2.0_dp) error stop "unit half weight c"
! Overflowing scaled weights return status 7 and do not write.
! Finite products, including those just below huge, are correctly rounded.
block
    integer(int64) :: bits, edge, gotb
    integer :: step, iscal
    real(dp) :: wv, expect
    xs(1) = 0.25_dp
    bits = transfer(1.0_dp, 0_int64)
    do step = -16, 260
        wv = transfer(bits + int(step, int64), 0.0_dp)
        iscal = 1024
        expect = scale(wv, iscal)
        ws(1) = wv
        xs(1) = 0.25_dp
        xkeep = xs(1)
        wkeep = ws(1)
        stc = gjp_scale_f(1_c_int, 1023.0_c_double, 0.0_c_double, &
                          -2.0_c_double, 2.0_c_double, xs, ws)
        if (.not. ieee_is_finite(expect)) then
            if (stc /= 7_c_int) error stop "pow2 1024 overflow status"
            if (xs(1) /= xkeep .or. ws(1) /= wkeep) error stop "pow2 1024 wrote"
        else
            if (stc /= 0_c_int .or. ws(1) /= expect) error stop "pow2 1024 finite"
        end if
    end do
    bits = transfer(2.0_dp, 0_int64)
    do step = -8, 16
        wv = transfer(bits + int(step, int64), 0.0_dp)
        iscal = 1023
        expect = scale(wv, iscal)
        ws(1) = wv
        xs(1) = 0.25_dp
        xkeep = xs(1)
        wkeep = ws(1)
        stc = gjp_scale_f(1_c_int, 1022.0_c_double, 0.0_c_double, &
                          -2.0_c_double, 2.0_c_double, xs, ws)
        if (.not. ieee_is_finite(expect)) then
            if (stc /= 7_c_int) error stop "pow2 1023 overflow status"
            if (xs(1) /= xkeep .or. ws(1) /= wkeep) error stop "pow2 1023 wrote"
        else
            if (stc /= 0_c_int .or. ws(1) /= expect) error stop "pow2 1023 finite"
        end if
    end do
    ws(1) = 1.0_dp
    xs(1) = 0.25_dp
    xkeep = xs(1)
    wkeep = ws(1)
    stc = gjp_scale_f(1_c_int, 511.0_c_double, 0.0_c_double, &
                      -4.0_c_double, 4.0_c_double, xs, ws)
    if (stc /= 7_c_int .or. xs(1) /= xkeep .or. ws(1) /= wkeep) &
        error stop "half 4 alpha 511 overflow"
    ws(1) = -1.0_dp
    xs(1) = 0.25_dp
    xkeep = xs(1)
    wkeep = ws(1)
    stc = gjp_scale_f(1_c_int, 1023.0_c_double, 0.0_c_double, &
                      -2.0_c_double, 2.0_c_double, xs, ws)
    if (stc /= 7_c_int .or. xs(1) /= xkeep .or. ws(1) /= wkeep) &
        error stop "negative weight overflow"
    ws(1) = nearest(1.0_dp, -1.0_dp)
    xs(1) = 0.25_dp
    call gauss_jacobi_scale(1, 1023.0_dp, 0.0_dp, -2.0_dp, 2.0_dp, xs, ws)
    if (ws(1) /= huge(1.0_dp)) error stop "predecessor of 1 scales to huge"
    wv = transfer(ior(ishft(int(z'3FEED709', int64), 32), &
                      int(z'438659DC', int64)), 0.0_dp)
    ws(1) = wv
    xs(1) = 0.25_dp
    call gauss_jacobi_scale(1, 3180.0_dp, 0.0_dp, 0.0_dp, 2.5_dp, xs, ws)
    gotb = transfer(ws(1), 0_int64)
    if (gotb /= ior(ishft(int(z'7FEFFFFF', int64), 32), int(z'FFFFFFF6', int64))) &
        error stop "alpha 3180 rounded weight"
    edge = ior(ishft(int(z'3FEED709', int64), 32), int(z'438659E5', int64))
    do step = -40, 16
        wv = transfer(edge + int(step, int64), 0.0_dp)
        ws(1) = wv
        xs(1) = 0.25_dp
        xkeep = xs(1)
        wkeep = ws(1)
        stc = gjp_scale_f(1_c_int, 3180.0_c_double, 0.0_c_double, &
                          0.0_c_double, 2.5_c_double, xs, ws)
        if (step < 0) then
            if (stc /= 0_c_int .or. .not. ieee_is_finite(ws(1))) &
                error stop "alpha 3180 finite side"
        else
            if (stc /= 7_c_int .or. xs(1) /= xkeep .or. ws(1) /= wkeep) &
                error stop "alpha 3180 overflow side"
        end if
    end do
    edge = ior(ishft(int(z'3FC702AE', int64), 32), int(z'4D1FB5D4', int64))
    do step = -201, 16
        wv = transfer(edge + int(step, int64), 0.0_dp)
        ws(1) = wv
        xs(1) = 0.25_dp
        xkeep = xs(1)
        wkeep = ws(1)
        stc = gjp_scale_f(1_c_int, 308.0_c_double, 0.0_c_double, &
                          0.0_c_double, 20.0_c_double, xs, ws)
        if (step < 0) then
            if (stc /= 0_c_int .or. .not. ieee_is_finite(ws(1))) &
                error stop "alpha 308 finite side"
        else
            if (stc /= 7_c_int .or. xs(1) /= xkeep .or. ws(1) /= wkeep) &
                error stop "alpha 308 overflow side"
        end if
    end do
    wv = transfer(ior(ishft(int(z'3FC702AE', int64), 32), &
                      int(z'4D1FB5B5', int64)), 0.0_dp)
    ws(1) = wv
    xs(1) = 0.25_dp
    stc = gjp_scale_f(1_c_int, 308.0_c_double, 0.0_c_double, &
                      0.0_c_double, 20.0_c_double, xs, ws)
    gotb = transfer(ws(1), 0_int64)
    if (stc /= 0_c_int) error stop "alpha 308 status"
    if (gotb /= ior(ishft(int(z'7FEFFFFF', int64), 32), int(z'FFFFFFD5', int64))) &
        error stop "alpha 308 rounded weight"
    wv = transfer(ior(ishft(int(z'3FEA8CE6', int64), 32), &
                      int(z'778183BB', int64)), 0.0_dp)
    ws(1) = wv
    xs(1) = -0.5_dp
    call gauss_jacobi_scale(1, 1750.0_dp, 0.0_dp, 0.0_dp, 3.0_dp, xs, ws)
    if (ws(1) /= huge(1.0_dp)) error stop "alpha 1750 rounds to huge"
end block
! Ties and signed zeros. Compare the bits, including the sign of zero.
block
    integer(int64) :: gotb, expb
    integer :: ie
    real(dp) :: wv
    xs(1) = 0.25_dp
    wv = transfer(ior(ishft(int(z'7F8948B0', int64), 32), int(z'FCD6E9E0', int64)), 0.0_dp)
    ws(1) = wv
    xkeep = xs(1)
    wkeep = ws(1)
    stc = gjp_scale_f(1_c_int, 0.0_c_double, 0.0_c_double, &
                      0.0_c_double, 162.0_c_double, xs, ws)
    if (stc /= 7_c_int .or. xs(1) /= xkeep .or. ws(1) /= wkeep) &
        error stop "half 81 overflow tie"
    ws(1) = -wv
    xs(1) = 0.25_dp
    xkeep = xs(1)
    wkeep = ws(1)
    stc = gjp_scale_f(1_c_int, 0.0_c_double, 0.0_c_double, &
                      0.0_c_double, 162.0_c_double, xs, ws)
    if (stc /= 7_c_int .or. xs(1) /= xkeep .or. ws(1) /= wkeep) &
        error stop "half 81 negative overflow tie"
    wv = transfer(ior(ishft(int(z'7FE99999', int64), 32), int(z'99999972', int64)), 0.0_dp)
    ws(1) = wv
    xs(1) = 0.25_dp
    stc = gjp_scale_f(1_c_int, 0.0_c_double, 0.0_c_double, &
                      0.0_c_double, 2.5_c_double, xs, ws)
    gotb = transfer(ws(1), 0_int64)
    expb = ior(ishft(int(z'7FEFFFFF', int64), 32), int(z'FFFFFFCE', int64))
    if (stc /= 0_c_int .or. gotb /= expb) error stop "half 1.25 interior tie"
    ws(1) = 1.0_dp
    xs(1) = 0.0_dp
    stc = gjp_scale_f(1_c_int, 106.5_c_double, 0.0_c_double, &
                      0.0_c_double, 2.0_c_double**(-9), xs, ws)
    if (stc /= 0_c_int .or. transfer(ws(1), 0_int64) /= 0_int64) &
        error stop "exact 2**-1075 is +0"
    ws(1) = -1.0_dp
    xs(1) = 0.0_dp
    stc = gjp_scale_f(1_c_int, 106.5_c_double, 0.0_c_double, &
                      0.0_c_double, 2.0_c_double**(-9), xs, ws)
    if (stc /= 0_c_int .or. transfer(ws(1), 0_int64) /= ishft(1_int64, 63)) &
        error stop "exact 2**-1075 is -0"
    ws(1) = 7.0_dp
    xs(1) = 0.0_dp
    stc = gjp_scale_f(1_c_int, 536.5_c_double, 0.0_c_double, &
                      0.0_c_double, 0.5_c_double, xs, ws)
    if (stc /= 0_c_int .or. transfer(ws(1), 0_int64) /= 4_int64) &
        error stop "7*2**-1075"
    ws(1) = -7.0_dp
    xs(1) = 0.0_dp
    stc = gjp_scale_f(1_c_int, 536.5_c_double, 0.0_c_double, &
                      0.0_c_double, 0.5_c_double, xs, ws)
    expb = ior(ishft(1_int64, 63), 4_int64)
    if (stc /= 0_c_int .or. transfer(ws(1), 0_int64) /= expb) &
        error stop "negative 7*2**-1075"
    ie = -1074
    ws(1) = 3.0_dp * scale(1.0_dp, ie)
    xs(1) = 0.0_dp
    stc = gjp_scale_f(1_c_int, 0.0_c_double, 0.0_c_double, &
                      0.0_c_double, 3.0_c_double, xs, ws)
    if (stc /= 0_c_int .or. transfer(ws(1), 0_int64) /= 4_int64) &
        error stop "1.5 subnormal tie"
    ws(1) = -3.0_dp * scale(1.0_dp, ie)
    xs(1) = 0.0_dp
    stc = gjp_scale_f(1_c_int, 0.0_c_double, 0.0_c_double, &
                      0.0_c_double, 3.0_c_double, xs, ws)
    if (stc /= 0_c_int .or. transfer(ws(1), 0_int64) /= expb) &
        error stop "negative 1.5 subnormal tie"
    ws(1) = 1.0_dp
    xs(1) = 0.0_dp
    stc = gjp_scale_f(1_c_int, -0.5_c_double, 0.0_c_double, &
                      0.0_c_double, nearest(0.0_c_double, 1.0_c_double), xs, ws)
    expb = ior(ishft(int(z'1E56A09E', int64), 32), int(z'667F3BCD', int64))
    if (stc /= 0_c_int .or. transfer(ws(1), 0_int64) /= expb) &
        error stop "flushed half square root"
    ws(1) = -1.0_dp
    xs(1) = 0.0_dp
    stc = gjp_scale_f(1_c_int, -0.5_c_double, 0.0_c_double, &
                      0.0_c_double, nearest(0.0_c_double, 1.0_c_double), xs, ws)
    expb = ior(ishft(int(z'9E56A09E', int64), 32), int(z'667F3BCD', int64))
    if (stc /= 0_c_int .or. transfer(ws(1), 0_int64) /= expb) &
        error stop "flushed half negative square root"
    wv = transfer(ior(ishft(int(z'3FB2C6D5', int64), 32), int(z'C2EFF05F', int64)), 0.0_dp)
    ws(1) = wv
    xs(1) = 0.0_dp
    stc = gjp_scale_f(1_c_int, 0.0_c_double, 0.0_c_double, &
                      0.0_c_double, 3.0_c_double, xs, ws)
    expb = ior(ishft(int(z'3FBC2A40', int64), 32), int(z'A467E88E', int64))
    if (stc /= 0_c_int .or. transfer(ws(1), 0_int64) /= expb) &
        error stop "Legendre [0,3] weight"
    wv = transfer(ior(ishft(int(z'3F9A0060', int64), 32), int(z'A8531FEC', int64)), 0.0_dp)
    ws(1) = wv
    xs(1) = 0.0_dp
    stc = gjp_scale_f(1_c_int, 0.0_c_double, 0.0_c_double, &
                      -2.0_c_double, 5.0_c_double, xs, ws)
    expb = ior(ishft(int(z'3FB6C054', int64), 32), int(z'9348BBEE', int64))
    if (stc /= 0_c_int .or. transfer(ws(1), 0_int64) /= expb) &
        error stop "Legendre [-2,5] weight"
    wv = transfer(ior(ishft(int(z'3F89465B', int64), 32), int(z'D31121F0', int64)), 0.0_dp)
    ws(1) = wv
    xs(1) = 0.0_dp
    stc = gjp_scale_f(1_c_int, 0.0_c_double, 0.0_c_double, &
                      -100.0_c_double, 100.0_c_double, xs, ws)
    expb = ior(ishft(int(z'3FF3BEF7', int64), 32), int(z'BCE56284', int64))
    if (stc /= 0_c_int .or. transfer(ws(1), 0_int64) /= expb) &
        error stop "Legendre [-100,100] weight"
    ! w * 3**2 on [0, 6] (alpha=1, p=2). Fraction ties-to-even: exact
    ! rational, remainder against half the denominator, even quotient kept.
    ! The product is the midpoint above even 0x1.0000000000006p+10.
    wv = transfer(ior(ishft(int(z'405C71C7', int64), 32), int(z'1C71C728', int64)), 0.0_dp)
    ws(1) = wv
    xs(1) = 0.25_dp
    stc = gjp_scale_f(1_c_int, 1.0_c_double, 0.0_c_double, &
                      0.0_c_double, 6.0_c_double, xs, ws)
    expb = ior(ishft(int(z'40900000', int64), 32), int(z'00000006', int64))
    if (stc /= 0_c_int .or. transfer(ws(1), 0_int64) /= expb) &
        error stop "half 3 square tie"
    ! w * 9**2 on [0, 18] is the overflow midpoint 2**1024 - 2**970.
    ! The same rounding is non-finite: status 7, both buffers unchanged.
    wv = transfer(ior(ishft(int(z'7F8948B0', int64), 32), int(z'FCD6E9E0', int64)), 0.0_dp)
    ws(1) = wv
    xs(1) = 0.25_dp
    xkeep = xs(1)
    wkeep = ws(1)
    stc = gjp_scale_f(1_c_int, 1.0_c_double, 0.0_c_double, &
                      0.0_c_double, 18.0_c_double, xs, ws)
    if (stc /= 7_c_int .or. xs(1) /= xkeep .or. ws(1) /= wkeep) &
        error stop "half 9 overflow tie"
end block
print '(A)', "scale endpoints and weight range OK"
end program test_scale_interval

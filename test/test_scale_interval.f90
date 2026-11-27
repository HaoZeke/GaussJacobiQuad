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
                                                          int(z'C00B0B8AD6203548', int64), int(z'40417FD9C31752FE', int64), &
                                                          int(z'BFEFFFFFFFFFFFFF', int64), int(z'C00B0B8AD6203543', int64), &
                                                          int(z'BFA95FA5F425D75E', int64), int(z'3FE5E952758AE35F', int64), &
                                                          int(z'BFEBB093A01F6698', int64), int(z'BF1B94A33920DC3A', int64), &
                                                          int(z'BF343DCC3E402080', int64), int(z'3FB00880ED342634', int64), &
                                                          int(z'BFEFFFFFFFFFFFFF', int64), int(z'BF343DCC3E402040', int64), &
                                                          int(z'7E37E43C8800759C', int64), int(z'7FEFFFFFFFFFFFFF', int64), &
                                                          int(z'BFEFFFFFFFFFFFFF', int64), int(z'7E37E43C8C00759C', int64), &
                                                          int(z'C1E0000000000000', int64), int(z'4540000000000000', int64), &
                                                          int(z'BFEFFFFFFFFFFFFF', int64), int(z'3E80000000000000', int64), &
                                                          int(z'0000000000000000', int64), int(z'4008000000000000', int64), &
                                                          int(z'BFEFFFFFFFFFFFFF', int64), int(z'3CA8000000000000', int64), &
                                                          int(z'BFF0000000000000', int64), int(z'4341C37937E08000', int64), &
                                                          int(z'BFEFFFFFFFFFFFFF', int64), int(z'BFDC790D903EFFFF', int64), &
                                                          int(z'4010322FE7216CB2', int64), int(z'4010322FE7216CB4', int64), &
                                                          int(z'3FD0000000000000', int64), int(z'4010322FE7216CB3', int64), &
                                                          int(z'C0F0000000000000', int64), int(z'4520000000000000', int64), &
                                                          int(z'BFEFFFFFFFFFFFFF', int64), int(z'41BFFF0000000000', int64), &
                                                          int(z'C11CBA6B387D1AE4', int64), int(z'412A81F5077331B5', int64), &
                                                          int(z'3FD160B51B47533C', int64), int(z'41173DC92735AF53', int64), &
                                                          int(z'C018002994F18826', int64), int(z'C018002994F18822', int64), &
                                                          int(z'3FD0000000000000', int64), int(z'C018002994F18824', int64), &
                                                          int(z'C240000000000000', int64), int(z'43D0000000000000', int64), &
                                                          int(z'BFEFFFFFFFFFFFFF', int64), int(z'C23FFFFFFF000000', int64), &
                                                          int(z'7FE1CCF385EBC8A0', int64), int(z'7FEFFFFFFFFFFFFF', int64), &
                                                          int(z'BFEFFFFFFFFFFFFF', int64), int(z'7FE1CCF385EBC8A0', int64), &
                                                          int(z'C12C80B97F64820E', int64), int(z'41237FE0AC94AC8B', int64), &
                                                          int(z'BFDEFB68CF907B80', int64), int(z'C1201ED9085CF85A', int64) &
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
print '(A)', "scale endpoints and weight range OK"
end program test_scale_interval

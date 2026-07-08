!> Tests for bogaert, glr, auto dispatch, and bogaert non-Legendre policy.
program test_auto_glr_bogaert
use GaussJacobiQuad, only: gauss_jacobi_rule, select_method_auto
use gjp_types, only: dp
use gjp_common, only: jacobi_zeroeth_moment
implicit none

integer :: me, n, i
real(dp) :: alpha, beta, tol, max_dx, max_dw, mu0, wsum
real(dp), allocatable :: x(:), w(:), xref(:), wref(:)
character(len=:), allocatable :: mauto, m
logical :: ok

me = this_image()
ok = .true.
tol = 1.0e-8_dp

! --- Auto policy strings ---
mauto = select_method_auto(64, 0.0_dp, 0.0_dp)
if (mauto /= "bogaert") then
    if (me == 1) print *, "FAIL auto Legendre large expected bogaert got ", mauto
    ok = .false.
end if
mauto = select_method_auto(16, 0.5_dp, 12.0_dp)
if (mauto /= "algo665") then
    if (me == 1) print *, "FAIL auto high-beta small-n expected algo665 got ", mauto
    ok = .false.
end if
mauto = select_method_auto(200, 0.5_dp, 0.5_dp)
if (mauto /= "rec") then
    if (me == 1) print *, "FAIL auto large mild expected rec got ", mauto
    ok = .false.
end if

! --- Bogaert vs algo665 (Legendre n=32) ---
n = 32
allocate (x(n), w(n), xref(n), wref(n))
call gauss_jacobi_rule(n, 0.0_dp, 0.0_dp, xref, wref, "algo665")
call gauss_jacobi_rule(n, 0.0_dp, 0.0_dp, x, w, "bogaert")
max_dx = maxval(abs(x - xref))
max_dw = maxval(abs(w - wref))
mu0 = jacobi_zeroeth_moment(0.0_dp, 0.0_dp)
wsum = sum(w)
if (me == 1) print '(A,ES12.5,A,ES12.5,A,ES12.5)', &
    "bogaert vs algo665 max|dx|=", max_dx, " max|dw|=", max_dw, " sum(w)-mu0=", wsum - mu0
if (max_dx > 1.0e-12_dp .or. max_dw > 1.0e-11_dp) ok = .false.  ! asymptotic tol
if (abs(wsum - mu0) > 1.0e-8_dp) ok = .false.
if (any([( .not. finite_dp(x(i)), i=1, n )]) .or. any([( .not. finite_dp(w(i)), i=1, n )])) ok = .false.
if (any(w <= 0.0_dp)) ok = .false.
deallocate (x, w, xref, wref)

! --- GLR vs algo665 (Jacobi) ---
n = 40
alpha = 0.5_dp
beta = 0.5_dp
allocate (x(n), w(n), xref(n), wref(n))
call gauss_jacobi_rule(n, alpha, beta, xref, wref, "algo665")
call gauss_jacobi_rule(n, alpha, beta, x, w, "glr")
max_dx = maxval(abs(x - xref))
max_dw = maxval(abs(w - wref))
if (me == 1) print '(A,ES12.5,A,ES12.5)', "glr vs algo665 max|dx|=", max_dx, " max|dw|=", max_dw
if (max_dx > tol .or. max_dw > 1.0e-6_dp) ok = .false.
deallocate (x, w, xref, wref)

! --- GLR vs algo665 (Legendre; phase-march path) ---
n = 64
allocate (x(n), w(n), xref(n), wref(n))
call gauss_jacobi_rule(n, 0.0_dp, 0.0_dp, xref, wref, "algo665")
call gauss_jacobi_rule(n, 0.0_dp, 0.0_dp, x, w, "glr")
max_dx = maxval(abs(x - xref))
max_dw = maxval(abs(w - wref))
if (me == 1) print '(A,ES12.5,A,ES12.5)', "glr Legendre vs algo665 max|dx|=", max_dx, " max|dw|=", max_dw
if (max_dx > tol .or. max_dw > 1.0e-6_dp) ok = .false.
! nodes must be strictly increasing (phase-chain order)
if (any([(x(i) >= x(i+1), i=1, n-1)])) then
    if (me == 1) print *, "FAIL glr nodes not strictly increasing"
    ok = .false.
end if
deallocate (x, w, xref, wref)


! --- Auto entry (omit method) two regimes ---
n = 48
allocate (x(n), w(n), xref(n), wref(n))
call gauss_jacobi_rule(n, 0.0_dp, 0.0_dp, x, w)  ! auto -> bogaert
call gauss_jacobi_rule(n, 0.0_dp, 0.0_dp, xref, wref, "bogaert")
max_dx = maxval(abs(x - xref))
if (me == 1) print '(A,ES12.5,A,A)', "auto Legendre vs forced bogaert max|dx|=", max_dx, &
    " method=", select_method_auto(n, 0.0_dp, 0.0_dp)
if (max_dx > 0.0_dp) ok = .false.  ! same path
deallocate (x, w, xref, wref)

n = 24
alpha = 1.0_dp
beta = 10.0_dp
allocate (x(n), w(n), xref(n), wref(n))
m = select_method_auto(n, alpha, beta)
call gauss_jacobi_rule(n, alpha, beta, x, w)  ! auto
call gauss_jacobi_rule(n, alpha, beta, xref, wref, m)
max_dx = maxval(abs(x - xref))
if (me == 1) print '(A,A,A,ES12.5)', "auto high-beta method=", m, " vs forced max|dx|=", max_dx
if (max_dx > 0.0_dp) ok = .false.
deallocate (x, w, xref, wref)

! --- Bogaert non-Legendre policy: must error ---
! exercised via a separate subprocess in CI; here we only document structure
if (me == 1) print *, "bogaert non-Legendre: error-stop policy (see bogaert_policy test harness)"

if (.not. ok) error stop "test_auto_glr_bogaert FAILED"
if (me == 1) print *, "test_auto_glr_bogaert PASSED"
sync all

contains
logical function finite_dp(z)
    real(dp), intent(in) :: z
    finite_dp = (z == z) .and. (abs(z) < huge(z))
end function finite_dp
end program test_auto_glr_bogaert

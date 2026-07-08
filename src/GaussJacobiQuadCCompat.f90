!> C-compatible bindings around the single-image method set.
!>
!> Primary entry: =gauss_jacobi_rule_c= (status codes; empty/"auto" → policy).
!> Status: 0 ok; 1 bad npts; 2 bad alpha; 3 bad beta; 4 unknown method;
!>         5 bogaert non-Legendre; 6 bogaert n < 21.
!> Does not =error stop= on policy failures so C/Python can handle them.
module GaussJacobiQuadCCompat
use iso_c_binding, only: c_double, c_int, c_char, c_null_char
use gjp_rec, only: gauss_jacobi_rec
use gjp_gw, only: gauss_jacobi_gw
use gjp_algo665, only: gauss_jacobi_algo665, gauss_jacobi_algo665_dc
use gjp_sturm, only: gauss_jacobi_sturm
use gjp_glr, only: gauss_jacobi_glr
use gjp_bogaert, only: gauss_jacobi_bogaert, bogaert_min_n
use gjp_auto, only: select_method_auto
use gjp_types, only: dp
implicit none

private
public :: gauss_jacobi_rule_c
public :: gauss_jacobi_rec_c, gauss_jacobi_gw_c, gauss_jacobi_algo665_c

integer(c_int), parameter, public :: GJP_OK = 0_c_int
integer(c_int), parameter, public :: GJP_ERR_NPTS = 1_c_int
integer(c_int), parameter, public :: GJP_ERR_ALPHA = 2_c_int
integer(c_int), parameter, public :: GJP_ERR_BETA = 3_c_int
integer(c_int), parameter, public :: GJP_ERR_METHOD = 4_c_int
integer(c_int), parameter, public :: GJP_ERR_BOGAERT_AB = 5_c_int
integer(c_int), parameter, public :: GJP_ERR_BOGAERT_N = 6_c_int

contains

!> Single C entry. =method_c= is a NUL-terminated C string; empty or "auto" selects policy.
function gauss_jacobi_rule_c(npts, alpha, beta, x, wts, method_c) result(status) &
    bind(C, name="gjp_rule_f")
    integer(c_int), intent(in), value :: npts
    real(c_double), intent(in), value :: alpha, beta
    real(c_double), intent(out) :: x(npts), wts(npts)
    character(kind=c_char), intent(in) :: method_c(*)
    integer(c_int) :: status
    character(len=:), allocatable :: m
    integer :: i, nlen

    status = GJP_OK
    if (npts <= 0) then
        status = GJP_ERR_NPTS
        return
    end if
    if (alpha <= -1.0_c_double) then
        status = GJP_ERR_ALPHA
        return
    end if
    if (beta <= -1.0_c_double) then
        status = GJP_ERR_BETA
        return
    end if

    nlen = 0
    do i = 1, 64
        if (method_c(i) == c_null_char) exit
        nlen = nlen + 1
    end do
    if (nlen == 0) then
        m = select_method_auto(int(npts), real(alpha, dp), real(beta, dp))
    else
        allocate (character(len=nlen) :: m)
        do i = 1, nlen
            m(i:i) = method_c(i)
        end do
        if (trim(m) == "auto" .or. len_trim(m) == 0) then
            deallocate (m)
            m = select_method_auto(int(npts), real(alpha, dp), real(beta, dp))
        else
            m = trim(m)
        end if
    end if

    status = dispatch_method(int(npts), real(alpha, dp), real(beta, dp), x, wts, m)
end function gauss_jacobi_rule_c

integer(c_int) function dispatch_method(npts, alpha, beta, x, wts, method) result(status)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(c_double), intent(out) :: x(npts), wts(npts)
    character(len=*), intent(in) :: method
    real(dp) :: xd(npts), wd(npts)

    status = GJP_OK
    select case (trim(method))
    case ("rec")
        call gauss_jacobi_rec(npts, alpha, beta, xd, wd)
    case ("gw")
        call gauss_jacobi_gw(npts, alpha, beta, xd, wd)
    case ("algo665")
        call gauss_jacobi_algo665(npts, alpha, beta, xd, wd)
    case ("algo665_dc")
        call gauss_jacobi_algo665_dc(npts, alpha, beta, xd, wd, use_caf=.false.)
    case ("sturm")
        call gauss_jacobi_sturm(npts, alpha, beta, xd, wd)
    case ("glr")
        call gauss_jacobi_glr(npts, alpha, beta, xd, wd)
    case ("bogaert")
        if (abs(alpha) > 0.0_dp .or. abs(beta) > 0.0_dp) then
            status = GJP_ERR_BOGAERT_AB
            return
        end if
        if (npts < bogaert_min_n) then
            status = GJP_ERR_BOGAERT_N
            return
        end if
        call gauss_jacobi_bogaert(npts, alpha, beta, xd, wd)
    case default
        status = GJP_ERR_METHOD
        return
    end select
    x = real(xd, c_double)
    wts = real(wd, c_double)
end function dispatch_method

subroutine gauss_jacobi_rec_c(npts, alpha, beta, x, wts) bind(C, name="gauss_jacobi_rec_c")
    integer(c_int), intent(in) :: npts
    real(c_double), intent(in) :: alpha, beta
    real(c_double), intent(out) :: x(npts), wts(npts)
    real(dp) :: xd(npts), wd(npts)
    call gauss_jacobi_rec(int(npts), real(alpha, dp), real(beta, dp), xd, wd)
    x = real(xd, c_double)
    wts = real(wd, c_double)
end subroutine gauss_jacobi_rec_c

subroutine gauss_jacobi_gw_c(npts, alpha, beta, x, wts) bind(C, name="gauss_jacobi_gw_c")
    integer(c_int), intent(in) :: npts
    real(c_double), intent(in) :: alpha, beta
    real(c_double), intent(out) :: x(npts), wts(npts)
    real(dp) :: xd(npts), wd(npts)
    call gauss_jacobi_gw(int(npts), real(alpha, dp), real(beta, dp), xd, wd)
    x = real(xd, c_double)
    wts = real(wd, c_double)
end subroutine gauss_jacobi_gw_c

subroutine gauss_jacobi_algo665_c(npts, alpha, beta, x, wts) bind(C, name="gauss_jacobi_algo665_c")
    integer(c_int), intent(in) :: npts
    real(c_double), intent(in) :: alpha, beta
    real(c_double), intent(out) :: x(npts), wts(npts)
    real(dp) :: xd(npts), wd(npts)
    call gauss_jacobi_algo665(int(npts), real(alpha, dp), real(beta, dp), xd, wd)
    x = real(xd, c_double)
    wts = real(wd, c_double)
end subroutine gauss_jacobi_algo665_c

end module GaussJacobiQuadCCompat

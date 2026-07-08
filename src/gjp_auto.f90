!> Method auto-selection policy for GaussJacobiQuad.
!>
!> Policy (documented, safe, not oracle-optimal):
!> 1. alpha=beta=0 and n >= 21 -> bogaert (iteration-free Legendre asymptotics)
!> 2. n >= 128 and max(|alpha|,|beta|) <= 5 -> rec (HT-style Newton; large n)
!> 3. n >= 64 and max(|alpha|,|beta|) <= 2 -> glr (Prufer + Newton)
!> 4. max(|alpha|,|beta|) > 8 or n <= 32 -> algo665 (robust ACM 655 / imtqlx)
!> 5. else -> sturm (GW via bisection; good general default)
module gjp_auto
use gjp_types, only: dp
use gjp_bogaert, only: bogaert_min_n
implicit none

private
public :: select_method_auto, method_is_legendre_only

contains

!> Return a deferred-length method string for auto policy.
function select_method_auto(npts, alpha, beta) result(method)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    character(len=:), allocatable :: method
    real(dp) :: abmax

    abmax = max(abs(alpha), abs(beta))

    if (abs(alpha) == 0.0_dp .and. abs(beta) == 0.0_dp .and. npts >= bogaert_min_n) then
        method = "bogaert"
    else if (npts >= 128 .and. abmax <= 5.0_dp) then
        method = "rec"
    else if (npts >= 64 .and. abmax <= 2.0_dp) then
        method = "glr"
    else if (abmax > 8.0_dp .or. npts <= 32) then
        method = "algo665"
    else
        method = "sturm"
    end if
end function select_method_auto

logical function method_is_legendre_only(method) result(leg)
    character(len=*), intent(in) :: method
    leg = (trim(method) == "bogaert")
end function method_is_legendre_only

end module gjp_auto

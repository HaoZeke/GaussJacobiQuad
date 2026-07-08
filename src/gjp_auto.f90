!> Method auto-selection policy for GaussJacobiQuad.
!>
!> Policy (documented, safe, not oracle-optimal). When num_images() > 1 the
!> CAF-capable names are preferred so multi-image runs actually partition work.
!>
!> Single-image:
!>  1. alpha=beta=0 and n >= 21 -> bogaert
!>  2. n >= 128 and max(|alpha|,|beta|) <= 5 -> rec
!>  3. n >= 64 and max(|alpha|,|beta|) <= 2 -> glr
!>  4. max(|alpha|,|beta|) > 8 or n <= 32 -> algo665
!>  5. else -> sturm
!>
!> Multi-image: same regimes with *_caf where available
!>  (bogaert_caf, rec_caf, glr_caf, sturm_caf; algo665 -> algo665_dc_caf for large n).
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
    logical :: multi

    abmax = max(abs(alpha), abs(beta))
    multi = (num_images() > 1)

    if (abs(alpha) == 0.0_dp .and. abs(beta) == 0.0_dp .and. npts >= bogaert_min_n) then
        if (multi) then
            method = "bogaert_caf"
        else
            method = "bogaert"
        end if
    else if (npts >= 128 .and. abmax <= 5.0_dp) then
        if (multi) then
            method = "rec_caf"
        else
            method = "rec"
        end if
    else if (npts >= 64 .and. abmax <= 2.0_dp) then
        if (multi) then
            method = "glr_caf"
        else
            method = "glr"
        end if
    else if (abmax > 8.0_dp .or. npts <= 32) then
        if (multi .and. npts >= 64) then
            method = "algo665_dc_caf"
        else if (multi) then
            method = "sturm_caf"
        else
            method = "algo665"
        end if
    else
        if (multi) then
            method = "sturm_caf"
        else
            method = "sturm"
        end if
    end if
end function select_method_auto

logical function method_is_legendre_only(method) result(leg)
    character(len=*), intent(in) :: method
    leg = (trim(method) == "bogaert" .or. trim(method) == "bogaert_caf")
end function method_is_legendre_only

end module gjp_auto

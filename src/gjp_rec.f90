! BEGIN_HEADER
! -----------------------------------------------------------------------------
! Gauss-Jacobi Quadrature Implementation
! Authors: Rohit Goswami <rgoswami[at]ieee.org>
! Source: GaussJacobiQuad Library
! License: MIT
! GitHub Repository: https://github.com/HaoZeke/GaussJacobiQuad
! Date: 2023-08-28
! Commit: c442f77
! -----------------------------------------------------------------------------
! This code is part of the GaussJacobiQuad library, providing an efficient
! implementation for Gauss-Jacobi quadrature nodes and weights computation.
! -----------------------------------------------------------------------------
! To cite this software:
! Rohit Goswami (2023). HaoZeke/GaussJacobiQuad: v0.1.0.
! Zenodo: https://doi.org/10.5281/ZENODO.8285112
! ---------------------------------------------------------------------
! END_HEADER
!> @brief This is derived from the recurrence relations in chebfun
!
!! @details This module implements the Gauss-Jacobi quadrature method for numerical integration.
!! It provides subroutines for computing the Jacobi polynomials, recurrence relations,
!! and evaluation of the polynomials and their derivatives.
!! The implementation is based on the recurrence method in chebfun (https://chebfun.org)
!! and the paper:
!! Hale and Townsend, Fast and accurate computation of Gauss–Legendre and
!! Gauss–Jacobi quadrature nodes and weights, SIAM J. Sci. Comp. 2013
!! DOI: 10.1137/120889873; docs/refs.bib key haleFastAccurateComputation2013
!!
!! Coarray Fortran: `gauss_jacobi_rec_caf` partitions independent Newton node
!! indices across images (`this_image` / `num_images`), then gathers via coarrays.
!! With one image the result matches the serial `gauss_jacobi_rec` path.
module gjp_rec
use gjp_types, only: dp
use gjp_constants, only: pi
implicit none

private
public :: gauss_jacobi_rec, gauss_jacobi_rec_caf
public :: caf_root_owner, caf_owns_index

contains

! This returns unsorted roots and weights (serial)
subroutine gauss_jacobi_rec(npts, alpha, beta, x, wts)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    real(dp), dimension(ceiling(npts/2._dp)) :: x1, ders1
    real(dp), dimension(npts/2) :: x2, ders2
    real(dp) :: ders(npts), C
    integer :: idx
    call recurrence(npts, ceiling(npts / 2._dp), alpha, beta, x1, ders1)
    call recurrence(npts, npts / 2, beta, alpha, x2, ders2)
    do idx = 1, npts / 2
        x(idx) = -x2(npts / 2 - idx + 1)
        ders(idx) = ders2(npts / 2 - idx + 1)
    end do
    do idx = 1, ceiling(npts / 2._dp)
        x(npts / 2 + idx) = x1(idx)
        ders(npts / 2 + idx) = ders1(idx)
    end do
    wts = 1.0_dp / ((1.0_dp - x**2) * ders**2)
    C = 2**(alpha + beta + 1) * exp(log_gamma(npts + alpha + 1) - &
                                    log_gamma(npts + alpha + beta + 1) + &
                                    log_gamma(npts + beta + 1) - log_gamma(npts + 1._dp))
    wts = wts * C
end subroutine gauss_jacobi_rec

!> Coarray-parallel recurrence path. Partition independent half-interval Newton
!> roots across images, gather with coarrays, then form nodes/weights.
!> One image owns every index (same math as serial; exercised under -fcoarray=single).
subroutine gauss_jacobi_rec_caf(npts, alpha, beta, x, wts)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    real(dp), dimension(ceiling(npts/2._dp)) :: x1, ders1
    real(dp), dimension(npts/2) :: x2, ders2
    real(dp) :: ders(npts), C
    integer :: idx

    call recurrence_caf(npts, ceiling(npts / 2._dp), alpha, beta, x1, ders1)
    call recurrence_caf(npts, npts / 2, beta, alpha, x2, ders2)
    do idx = 1, npts / 2
        x(idx) = -x2(npts / 2 - idx + 1)
        ders(idx) = ders2(npts / 2 - idx + 1)
    end do
    do idx = 1, ceiling(npts / 2._dp)
        x(npts / 2 + idx) = x1(idx)
        ders(npts / 2 + idx) = ders1(idx)
    end do
    wts = 1.0_dp / ((1.0_dp - x**2) * ders**2)
    C = 2**(alpha + beta + 1) * exp(log_gamma(npts + alpha + 1) - &
                                    log_gamma(npts + alpha + beta + 1) + &
                                    log_gamma(npts + beta + 1) - log_gamma(npts + 1._dp))
    wts = wts * C
end subroutine gauss_jacobi_rec_caf

!> Image that owns 1-based index `i` under round-robin striding (1..nimg).
pure integer function caf_root_owner(i, nimg) result(owner)
    integer, intent(in) :: i, nimg
    if (nimg <= 1) then
        owner = 1
    else
        owner = mod(i - 1, nimg) + 1
    end if
end function caf_root_owner

!> True if image `me` owns index `i`.
pure logical function caf_owns_index(i, me, nimg) result(owns)
    integer, intent(in) :: i, me, nimg
    owns = (caf_root_owner(i, nimg) == me)
end function caf_owns_index

subroutine recurrence(npts, n2, alpha, beta, x, PP)
    integer, intent(in) :: npts, n2
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(n2), PP(n2)
    real(dp) :: dx(n2), P(n2)
    integer :: r(n2), l, i
    real(dp) :: C, T

    if (n2 <= 0) return

    do i = 1, n2
        r(i) = n2 - i + 1
    end do

    do i = 1, n2
        C = (2 * r(i) + alpha - 0.5_dp) * pi / (2 * npts + alpha + beta + 1)
     T = C + 1 / (2 * npts + alpha + beta + 1)**2 * ((0.25_dp - alpha**2) / tan(0.5_dp * C) - (0.25_dp - beta**2) * tan(0.5_dp * C))
        x(i) = cos(T)
    end do

    dx = 1.0_dp
    l = 0

    do while (maxval(abs(dx)) > sqrt(epsilon(1.0_dp)) / 1000 .and. l < 10)
        l = l + 1
        call eval_jacobi_poly(x, npts, alpha, beta, P, PP)
        dx = -P / PP
        x = x + dx
    end do

    call eval_jacobi_poly(x, npts, alpha, beta, P, PP)
end subroutine

!> Partition Newton roots across images; gather with coarrays after sync.
subroutine recurrence_caf(npts, n2, alpha, beta, x, PP)
    integer, intent(in) :: npts, n2
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(n2), PP(n2)
    real(dp), allocatable :: x_caf(:) [:], PP_caf(:) [:]
    integer :: me, nimg, i, img, r_i, l
    real(dp) :: C, T, dx, x_i, P_val, PP_val
    real(dp) :: xv(1), Pv(1), PPv(1)

    if (n2 <= 0) return

    me = this_image()
    nimg = num_images()

    ! Always use coarrays + image ownership so the CAF path is exercised even
    ! under -fcoarray=single (nimg==1 owns every index). Multi-image partitions
    ! indices round-robin; gather after sync all.
    allocate (x_caf(n2) [*], PP_caf(n2) [*])
    x_caf = 0.0_dp
    PP_caf = 0.0_dp

    do i = 1, n2
        if (.not. caf_owns_index(i, me, nimg)) cycle

        r_i = n2 - i + 1
        C = (2 * r_i + alpha - 0.5_dp) * pi / (2 * npts + alpha + beta + 1)
        T = C + 1 / (2 * npts + alpha + beta + 1)**2 * &
            ((0.25_dp - alpha**2) / tan(0.5_dp * C) - (0.25_dp - beta**2) * tan(0.5_dp * C))
        x_i = cos(T)

        dx = 1.0_dp
        l = 0
        do while (abs(dx) > sqrt(epsilon(1.0_dp)) / 1000 .and. l < 10)
            l = l + 1
            xv(1) = x_i
            call eval_jacobi_poly(xv, npts, alpha, beta, Pv, PPv)
            P_val = Pv(1)
            PP_val = PPv(1)
            dx = -P_val / PP_val
            x_i = x_i + dx
        end do
        xv(1) = x_i
        call eval_jacobi_poly(xv, npts, alpha, beta, Pv, PPv)
        x_caf(i) = x_i
        PP_caf(i) = PPv(1)
    end do

    sync all

    ! Image 1 gathers; full vectors redistributed from image 1
    if (me == 1) then
        do i = 1, n2
            img = caf_root_owner(i, nimg)
            if (img == 1) then
                x(i) = x_caf(i)
                PP(i) = PP_caf(i)
            else
                x(i) = x_caf(i) [img]
                PP(i) = PP_caf(i) [img]
            end if
        end do
        x_caf = x
        PP_caf = PP
    end if
    sync all
    x = x_caf(:) [1]
    PP = PP_caf(:) [1]

    sync all
    deallocate (x_caf, PP_caf)
end subroutine recurrence_caf

subroutine eval_jacobi_poly(x, npts, alpha, beta, P, Pp)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(in) :: x(:)
    real(dp), dimension(size(x)), intent(out) :: P, Pp
    real(dp), dimension(size(x)) :: Pm1, Ppm1, Pa1, Ppa1
    integer :: k, i
    real(dp) :: A_val, B_val, C_val, D_val

    P = 0.5_dp * (alpha - beta + (alpha + beta + 2) * x)
    Pm1 = 1.0_dp
    Pp = 0.5_dp * (alpha + beta + 2)
    Ppm1 = 0.0_dp

    if (npts == 0) then
        P = Pm1
        Pp = Ppm1
    end if

    do k = 1, npts - 1
        A_val = 2 * (k + 1) * (k + alpha + beta + 1) * (2 * k + alpha + beta)
        B_val = (2 * k + alpha + beta + 1) * (alpha**2 - beta**2)
        C_val = product([(2 * k + alpha + beta + i, i=0, 2)])
        D_val = 2 * (k + alpha) * (k + beta) * (2 * k + alpha + beta + 2)

        Pa1 = ((B_val + C_val * x) * P - D_val * Pm1) / A_val
        Ppa1 = ((B_val + C_val * x) * Pp + C_val * P - D_val * Ppm1) / A_val

        Pm1 = P
        P = Pa1
        Ppm1 = Pp
        Pp = Ppa1
    end do
end subroutine

end module gjp_rec

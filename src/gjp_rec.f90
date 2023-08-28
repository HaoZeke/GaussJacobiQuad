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
module gjp_rec
use gjp_types, only: dp
use gjp_constants, only: pi
implicit none
contains

! This returns unsorted roots and weights
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

subroutine recurrence(npts, n2, alpha, beta, x, PP)
    integer, intent(in) :: npts, n2
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(n2), PP(n2)
    real(dp) :: dx(n2), P(n2)
    integer :: r(n2), l, i
    real(dp) :: C, T

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

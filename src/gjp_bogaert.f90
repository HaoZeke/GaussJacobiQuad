!> Bogaert-style iteration-free Gauss-Legendre (alpha = beta = 0).
!>
!> Nodes from the asymptotic expansion of Bogaert (SIAM J. Sci. Comput. 2014)
!> about zeros of Bessel J0, matching the structure used in chebfun legpts
!> ('ASY'). Weights from the exact Gauss formula
!>   w = 2 / ((1-x^2) P_n'(x)^2)
!> with a three-term recurrence for P_n' at those nodes (iteration-free root
!> find; one recurrence pass for weights).
!>
!> Requires n >= 21 and alpha = beta = 0. Callers must reject/fallback otherwise.
!> DOI: 10.1137/140954969 (Bogaert 2014); see docs/refs.bib key bogaertIterationfreeComputationGausslegendre2014
module gjp_bogaert
use gjp_types, only: dp
use gjp_constants, only: pi
implicit none

private
public :: gauss_jacobi_bogaert, bogaert_min_n

integer, parameter :: bogaert_min_n = 21

contains

subroutine gauss_jacobi_bogaert(npts, alpha, beta, x, wts)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    integer :: m, i
    real(dp) :: a(npts), vn, vn2, vn4, vn6
    real(dp) :: u, u2, u4, ai, ai2, ai3, ai5, t
    real(dp) :: F0, F1, F2, F3, R30, R31, R32, R33, R35
    real(dp) :: p, pp, xi

    if (abs(alpha) > 0.0_dp .or. abs(beta) > 0.0_dp) then
        error stop "bogaert: requires alpha=beta=0 (Gauss-Legendre only)"
    end if
    if (npts < bogaert_min_n) then
        error stop "bogaert: requires n >= 21 for double-precision asymptotics"
    end if

    m = (npts + 1) / 2
    call bessel_j0_zeros(m, a)

    vn = 1.0_dp / (real(npts, dp) + 0.5_dp)
    do i = 1, m
        a(i) = a(i) * vn
    end do

    vn2 = vn * vn
    vn4 = vn2 * vn2
    vn6 = vn4 * vn2

    ! Bogaert/chebfun expansion: t = F0 + F1*vn^2 + F2*vn^4 + F3*vn^6
    ! then positive half-nodes x = cos(t); we store ascending full set.
    do i = 1, m
        ai = a(i)
        u = cos(ai) / sin(ai)  ! cot
        u2 = u * u
        u4 = u2 * u2
        ai2 = ai * ai
        ai3 = ai2 * ai
        ai5 = ai2 * ai3

        F0 = ai
        F1 = (u * ai - 1.0_dp) / (8.0_dp * ai)

        if (npts < 10000) then
            F2 = (6.0_dp * ai2 * (1.0_dp + u2) + 25.0_dp - u * (31.0_dp * u2 + 33.0_dp) * ai3) &
                 / (384.0_dp * ai3)
        else
            F2 = 0.0_dp
        end if

        if (npts < 1000) then
            R30 = u * (2595.0_dp + 6350.0_dp * u2 + 3779.0_dp * u4) / 15360.0_dp
            R31 = -(31.0_dp * u2 + 11.0_dp) / 1024.0_dp
            R32 = u / 512.0_dp
            R33 = -25.0_dp / 3072.0_dp
            R35 = -1073.0_dp / 5120.0_dp
            F3 = R30 + R35 / ai5 + (1.0_dp + u2) * (R31 / ai + R32 / ai2 + R33 / ai3)
        else
            F3 = 0.0_dp
        end if

        t = F0 + F1 * vn2 + F2 * vn4 + F3 * vn6
        ! half-array: largest positive first in chebfun; we fill left (negative) end
        x(i) = -cos(t)
    end do

    ! Weights via exact Gauss formula + Legendre recurrence at asymptotic nodes
    do i = 1, m
        xi = x(i)
        call legendre_p_pp(npts, xi, p, pp)
        wts(i) = 2.0_dp / ((1.0_dp - xi * xi) * pp * pp)
    end do

    do i = 1, m
        x(npts + 1 - i) = -x(i)
        wts(npts + 1 - i) = wts(i)
    end do
    if (mod(npts, 2) == 1) then
        x(m) = 0.0_dp
        ! recompute centre weight for odd n (x=0)
        call legendre_p_pp(npts, 0.0_dp, p, pp)
        wts(m) = 2.0_dp / (pp * pp)
    end if
end subroutine gauss_jacobi_bogaert

!> First m positive zeros of J0 via McMahon starter + Newton with J0/J1.
subroutine bessel_j0_zeros(m, z)
    integer, intent(in) :: m
    real(dp), intent(out) :: z(m)
    integer :: s, it
    real(dp) :: beta, j0, j1, dz

    do s = 1, m
        beta = (real(s, dp) - 0.25_dp) * pi
        z(s) = beta - 1.0_dp / (8.0_dp * beta) &
               - 7.0_dp / (384.0_dp * beta**3) &
               - 83.0_dp / (5120.0_dp * beta**5) &
               - 6949.0_dp / (1658880.0_dp * beta**7)
        ! Polish: J0'(z) = -J1(z)
        do it = 1, 8
            j0 = bessel_j0(z(s))
            j1 = bessel_j1(z(s))
            if (abs(j1) < epsilon(1.0_dp)) exit
            dz = j0 / j1
            z(s) = z(s) + dz
            if (abs(dz) < 1.0e-15_dp * z(s)) exit
        end do
    end do
end subroutine bessel_j0_zeros

!> P_n(x) and P_n'(x) for Legendre via Bonnet recurrence.
subroutine legendre_p_pp(n, x, p, pp)
    integer, intent(in) :: n
    real(dp), intent(in) :: x
    real(dp), intent(out) :: p, pp
    real(dp) :: pm1, ppm1, pk, ppk
    integer :: k

    if (n == 0) then
        p = 1.0_dp
        pp = 0.0_dp
        return
    end if
    pm1 = 1.0_dp
    p = x
    ppm1 = 0.0_dp
    pp = 1.0_dp
    do k = 1, n - 1
        pk = ((2 * k + 1) * x * p - k * pm1) / real(k + 1, dp)
        ppk = ((2 * k + 1) * (p + x * pp) - k * ppm1) / real(k + 1, dp)
        pm1 = p
        p = pk
        ppm1 = pp
        pp = ppk
    end do
end subroutine legendre_p_pp

end module gjp_bogaert

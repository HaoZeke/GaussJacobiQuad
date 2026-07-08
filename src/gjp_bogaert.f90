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
!> CAF: half-node indices are independent → partition across images.
!> DOI: 10.1137/140954969 (Bogaert 2014); see docs/refs.bib key bogaertIterationfreeComputationGausslegendre2014
module gjp_bogaert
use gjp_types, only: dp
use gjp_constants, only: pi
use gjp_rec, only: caf_root_owner, caf_owns_index
implicit none

private
public :: gauss_jacobi_bogaert, gauss_jacobi_bogaert_caf, bogaert_min_n

integer, parameter :: bogaert_min_n = 21

contains

subroutine gauss_jacobi_bogaert(npts, alpha, beta, x, wts)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    integer :: m, i
    real(dp) :: xi, wi

    call bogaert_check(npts, alpha, beta)
    m = (npts + 1) / 2
    do i = 1, m
        call bogaert_half_node(npts, i, xi, wi)
        x(i) = xi
        wts(i) = wi
        x(npts + 1 - i) = -xi
        wts(npts + 1 - i) = wi
    end do
    if (mod(npts, 2) == 1) then
        x(m) = 0.0_dp
        block
            real(dp) :: p0, pp0
            call legendre_p_pp(npts, 0.0_dp, p0, pp0)
            wts(m) = 2.0_dp / (pp0 * pp0)
        end block
    end if
end subroutine gauss_jacobi_bogaert

!> CAF: partition half-node indices (1..m) across images; gather full rule.
subroutine gauss_jacobi_bogaert_caf(npts, alpha, beta, x, wts)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    real(dp), allocatable :: x_c(:)[:], w_c(:)[:]
    integer :: me, nimg, m, i, img
    real(dp) :: xi, wi, p, pp

    call bogaert_check(npts, alpha, beta)
    me = this_image()
    nimg = num_images()
    m = (npts + 1) / 2

    allocate (x_c(npts)[*], w_c(npts)[*])
    x_c = 0.0_dp
    w_c = 0.0_dp

    do i = 1, m
        if (.not. caf_owns_index(i, me, nimg)) cycle
        call bogaert_half_node(npts, i, xi, wi)
        x_c(i) = xi
        w_c(i) = wi
        if (i /= npts + 1 - i) then
            x_c(npts + 1 - i) = -xi
            w_c(npts + 1 - i) = wi
        end if
    end do

    ! Odd n: centre node is half-index m with x=0; recompute weight if we own m
    if (mod(npts, 2) == 1 .and. caf_owns_index(m, me, nimg)) then
        x_c(m) = 0.0_dp
        call legendre_p_pp(npts, 0.0_dp, p, pp)
        w_c(m) = 2.0_dp / (pp * pp)
    end if

    sync all
    if (me == 1) then
        do i = 1, npts
            img = caf_root_owner(min(i, npts + 1 - i), nimg)
            ! Full index i may be the mirror of half-index j=min(i,n+1-i)
            if (img == 1) then
                x(i) = x_c(i)
                wts(i) = w_c(i)
            else
                x(i) = x_c(i)[img]
                wts(i) = w_c(i)[img]
            end if
        end do
        x_c = x
        w_c = wts
    end if
    sync all
    x = x_c(:)[1]
    wts = w_c(:)[1]
    sync all
    deallocate (x_c, w_c)
end subroutine gauss_jacobi_bogaert_caf

subroutine bogaert_check(npts, alpha, beta)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    if (abs(alpha) > 0.0_dp .or. abs(beta) > 0.0_dp) then
        error stop "bogaert: requires alpha=beta=0 (Gauss-Legendre only)"
    end if
    if (npts < bogaert_min_n) then
        error stop "bogaert: requires n >= 21 for double-precision asymptotics"
    end if
end subroutine bogaert_check

!> One left half-node (index i of m=(n+1)/2): asymptotic + weight at that node.
subroutine bogaert_half_node(npts, i, xi, wi)
    integer, intent(in) :: npts, i
    real(dp), intent(out) :: xi, wi
    real(dp) :: a, vn, vn2, vn4, vn6
    real(dp) :: u, u2, u4, ai, ai2, ai3, ai5, t
    real(dp) :: F0, F1, F2, F3, R30, R31, R32, R33, R35
    real(dp) :: p, pp

    call bessel_j0_zero_s(i, a)
    vn = 1.0_dp / (real(npts, dp) + 0.5_dp)
    a = a * vn
    vn2 = vn * vn
    vn4 = vn2 * vn2
    vn6 = vn4 * vn2

    ai = a
    u = cos(ai) / sin(ai)
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
    xi = -cos(t)
    call legendre_p_pp(npts, xi, p, pp)
    wi = 2.0_dp / ((1.0_dp - xi * xi) * pp * pp)
end subroutine bogaert_half_node

!> s-th positive zero of J0 (McMahon + Newton).
subroutine bessel_j0_zero_s(s, z)
    integer, intent(in) :: s
    real(dp), intent(out) :: z
    integer :: it
    real(dp) :: beta, j0, j1, dz

    beta = (real(s, dp) - 0.25_dp) * pi
    z = beta - 1.0_dp / (8.0_dp * beta) &
        - 7.0_dp / (384.0_dp * beta**3) &
        - 83.0_dp / (5120.0_dp * beta**5) &
        - 6949.0_dp / (1658880.0_dp * beta**7)
    do it = 1, 8
        j0 = bessel_j0(z)
        j1 = bessel_j1(z)
        if (abs(j1) < epsilon(1.0_dp)) exit
        dz = j0 / j1
        z = z + dz
        if (abs(dz) < 1.0e-15_dp * z) exit
    end do
end subroutine bessel_j0_zero_s

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

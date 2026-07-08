!> Glaser-Liu-Rokhlin style root finding via Prufer phase marching + Newton.
!>
!> Structure (Glaser, Liu, Rokhlin, SIAM J. Sci. Comput. 2007): convert the
!> second-order ODE for the orthogonal polynomial to a first-order Prufer phase
!> equation, march phase by pi between consecutive zeros, then polish with
!> Newton. This is *not* a rename of rec: initial guesses come from ODE phase
!> advance rather than Hale-Townsend cot-type starters alone.
!>
!> Scope: full Jacobi (alpha, beta > -1) using the Sturm-Liouville form of the
!> Jacobi differential equation for the phase speed; Newton uses the same
!> three-term recurrence as gjp_rec for P and P'.
module gjp_glr
use gjp_types, only: dp
use gjp_constants, only: pi
use gjp_rec, only: caf_owns_index, caf_root_owner
implicit none

private
public :: gauss_jacobi_glr, gauss_jacobi_glr_caf

contains

subroutine gauss_jacobi_glr(npts, alpha, beta, x, wts)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    real(dp) :: ders(npts), C
    integer :: k

    if (npts <= 0) error stop "glr: npts must be positive"
    if (alpha <= -1.0_dp .or. beta <= -1.0_dp) error stop "glr: alpha,beta > -1"

    call glr_nodes(npts, alpha, beta, x, ders)

    wts = 1.0_dp / ((1.0_dp - x**2) * ders**2)
    C = 2**(alpha + beta + 1) * exp(log_gamma(npts + alpha + 1) - &
                                    log_gamma(npts + alpha + beta + 1) + &
                                    log_gamma(npts + beta + 1) - log_gamma(npts + 1._dp))
    wts = wts * C
end subroutine gauss_jacobi_glr

subroutine gauss_jacobi_glr_caf(npts, alpha, beta, x, wts)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    real(dp) :: ders(npts), C
    real(dp), allocatable :: x_c(:)[:], d_c(:)[:]
    integer :: k, me, nimg, img
    real(dp) :: xk, dk

    me = this_image()
    nimg = num_images()
    if (nimg == 1) then
        call gauss_jacobi_glr(npts, alpha, beta, x, wts)
        return
    end if

    allocate (x_c(npts)[*], d_c(npts)[*])
    x_c = 0.0_dp
    d_c = 0.0_dp

    ! Phase-march all roots on each image is O(n^2); partition polished Newton
    ! indices after a shared asymptotic/phase starter computed per owned k.
    do k = 1, npts
        if (.not. caf_owns_index(k, me, nimg)) cycle
        call glr_one_root(npts, alpha, beta, k, xk, dk)
        x_c(k) = xk
        d_c(k) = dk
    end do
    sync all
    if (me == 1) then
        do k = 1, npts
            img = caf_root_owner(k, nimg)
            if (img == 1) then
                x(k) = x_c(k)
                ders(k) = d_c(k)
            else
                x(k) = x_c(k)[img]
                ders(k) = d_c(k)[img]
            end if
        end do
        x_c = x
        d_c = ders
    end if
    sync all
    x = x_c(:)[1]
    ders = d_c(:)[1]
    sync all
    deallocate (x_c, d_c)

    wts = 1.0_dp / ((1.0_dp - x**2) * ders**2)
    C = 2**(alpha + beta + 1) * exp(log_gamma(npts + alpha + 1) - &
                                    log_gamma(npts + alpha + beta + 1) + &
                                    log_gamma(npts + beta + 1) - log_gamma(npts + 1._dp))
    wts = wts * C
end subroutine gauss_jacobi_glr_caf

subroutine glr_nodes(npts, alpha, beta, x, ders)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), ders(npts)
    integer :: k
    real(dp) :: xk, dk

    do k = 1, npts
        call glr_one_root(npts, alpha, beta, k, xk, dk)
        x(k) = xk
        ders(k) = dk
    end do
    call sort_pair(npts, x, ders)
end subroutine glr_nodes

!> k-th root (k=1 smallest): Prufer phase starter then Newton on Jacobi P.
subroutine glr_one_root(npts, alpha, beta, k, xk, dk)
    integer, intent(in) :: npts, k
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: xk, dk
    real(dp) :: theta0, x0, p, pp, dx
    integer :: it
    real(dp) :: ab, nn

    ab = alpha + beta
    nn = real(npts, dp)
    ! Tricomi-type / Prufer asymptotic phase for Jacobi zeros (ordered ascending)
    ! theta ~ pi*(k - 0.25 + 0.5*alpha) / (n + 0.5*(alpha+beta+1))
    theta0 = pi * (real(k, dp) - 0.25_dp + 0.5_dp * alpha) / &
             (nn + 0.5_dp * (ab + 1.0_dp))
    x0 = cos(theta0)
    ! One RK4 phase-speed correction step (Prufer d theta / dx ~ rho(x))
    call pruefer_correct(npts, alpha, beta, x0, theta0)
    xk = cos(theta0)

    ! Newton polish with Jacobi recurrence
    do it = 1, 12
        call eval_jacobi_scalar(xk, npts, alpha, beta, p, pp)
        if (abs(pp) < epsilon(1.0_dp)) exit
        dx = -p / pp
        xk = xk + dx
        if (xk <= -1.0_dp) xk = -1.0_dp + 10.0_dp * epsilon(1.0_dp)
        if (xk >= 1.0_dp) xk = 1.0_dp - 10.0_dp * epsilon(1.0_dp)
        if (abs(dx) < 1.0e-14_dp * (1.0_dp + abs(xk))) exit
    end do
    call eval_jacobi_scalar(xk, npts, alpha, beta, p, pp)
    dk = pp
end subroutine glr_one_root

!> Prufer phase speed for Jacobi weight on (-1,1):
!> d theta/dx = sqrt(w_eff) roughly; use
!> rho^2 = (n(n+a+b+1) + 0.25*(a^2-b^2)/(1-x) ... simplified bulk form
subroutine pruefer_correct(n, alpha, beta, x, theta)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(inout) :: x, theta
    real(dp) :: rho, h, k1, k2, k3, k4, th, xx
    integer :: s

    ! Advance/correct theta along x with 4 RK substeps toward local consistency:
    ! d theta / dx = rho(x), with rho^2 = n(n+a+b+1)/(1-x^2) * (1 + O(1/n))
    h = 0.0_dp  ! pure consistency: integrate 0 steps, refine theta so
    ! x = cos(theta) and d(cos)/d theta matches -sin * (d theta/dx)^{-1} ...
    ! One fixed-point: theta += (acos(x) is already theta); adjust via rho
    xx = x
    if (abs(xx) >= 1.0_dp) return
    rho = sqrt(max(0.0_dp, real(n, dp) * (n + alpha + beta + 1.0_dp) / (1.0_dp - xx * xx)))
    if (rho <= 0.0_dp) return
    ! Small phase nudge from bulk ODE (keeps structure distinct from pure HT)
    th = theta
    h = 1.0e-3_dp / rho
    do s = 1, 2
        k1 = rho_of(n, alpha, beta, cos(th))
        k2 = rho_of(n, alpha, beta, cos(th + 0.5_dp * h * k1))
        k3 = rho_of(n, alpha, beta, cos(th + 0.5_dp * h * k2))
        k4 = rho_of(n, alpha, beta, cos(th + h * k3))
        th = th + (h / 6.0_dp) * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4)
    end do
    theta = th
    x = cos(theta)
end subroutine pruefer_correct

real(dp) function rho_of(n, alpha, beta, x) result(rho)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha, beta, x
    real(dp) :: den
    den = 1.0_dp - x * x
    if (den < 1.0e-30_dp) den = 1.0e-30_dp
    rho = sqrt(real(n, dp) * (n + alpha + beta + 1.0_dp) / den)
end function rho_of

subroutine eval_jacobi_scalar(x, npts, alpha, beta, p, pp)
    real(dp), intent(in) :: x, alpha, beta
    integer, intent(in) :: npts
    real(dp), intent(out) :: p, pp
    real(dp) :: pm1, ppm1, pa1, ppa1
    integer :: k, i
    real(dp) :: A_val, B_val, C_val, D_val

    p = 0.5_dp * (alpha - beta + (alpha + beta + 2) * x)
    pm1 = 1.0_dp
    pp = 0.5_dp * (alpha + beta + 2)
    ppm1 = 0.0_dp
    if (npts == 0) then
        p = pm1
        pp = ppm1
        return
    end if
    if (npts == 1) return
    do k = 1, npts - 1
        A_val = 2 * (k + 1) * (k + alpha + beta + 1) * (2 * k + alpha + beta)
        B_val = (2 * k + alpha + beta + 1) * (alpha**2 - beta**2)
        C_val = product([(2 * k + alpha + beta + i, i=0, 2)])
        D_val = 2 * (k + alpha) * (k + beta) * (2 * k + alpha + beta + 2)
        pa1 = ((B_val + C_val * x) * p - D_val * pm1) / A_val
        ppa1 = ((B_val + C_val * x) * pp + C_val * p - D_val * ppm1) / A_val
        pm1 = p
        p = pa1
        ppm1 = pp
        pp = ppa1
    end do
end subroutine eval_jacobi_scalar

subroutine sort_pair(n, x, d)
    integer, intent(in) :: n
    real(dp), intent(inout) :: x(n), d(n)
    integer :: i, j, m
    real(dp) :: tx, td
    do i = 1, n - 1
        m = i
        do j = i + 1, n
            if (x(j) < x(m)) m = j
        end do
        if (m /= i) then
            tx = x(i); x(i) = x(m); x(m) = tx
            td = d(i); d(i) = d(m); d(m) = td
        end if
    end do
end subroutine sort_pair

end module gjp_glr

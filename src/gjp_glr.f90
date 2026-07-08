!> Glaser–Liu–Rokhlin (GLR) root finding via Prüfer phase marching + Newton.
!>
!> Structure (Glaser, Liu, Rokhlin, SIAM J. Sci. Comput. 2007): convert the
!> second-order ODE for the orthogonal polynomial to a first-order Prüfer phase
!> equation, *march the phase by π between consecutive zeros*, then polish with
!> Newton. This is deliberately *not* independent asymptotic starters per index
!> (that is the rec / Hale–Townsend class): roots are chained — each guess is the
!> result of integrating dθ/dx = ρ(x) through one half-oscillation from the
!> previous polished zero.
!>
!> Scope: full Jacobi (α, β > −1). Local wave number ρ from the bulk Liouville–
!> Green / SL form of the Jacobi equation; Newton uses the three-term recurrence
!> for P and P′.
!>
!> Serial =glr=: true sequential Prüfer phase-by-π march (GLR).
!> Multi-image =glr_caf=: independent k-th root asymptotics + Newton (HT-style
!> partition of root indices). Phase marching is inherently serial; CAF uses the
!> same Newton polish / weight formula with parallel starters.
!> DOI: 10.1137/06067016x (Glaser–Liu–Rokhlin 2007); see docs/refs.bib key glaserFastAlgorithmCalculation2007
module gjp_glr
use gjp_types, only: dp
use gjp_constants, only: pi
use gjp_rec, only: caf_root_owner, caf_owns_index
implicit none

private
public :: gauss_jacobi_glr, gauss_jacobi_glr_caf

integer, parameter :: newton_max_it = 16
integer, parameter :: rk_steps_per_pi = 48
real(dp), parameter :: x_eps = 1.0e-14_dp

contains

subroutine gauss_jacobi_glr(npts, alpha, beta, x, wts)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    real(dp) :: ders(npts), C

    if (npts <= 0) error stop "glr: npts must be positive"
    if (alpha <= -1.0_dp .or. beta <= -1.0_dp) error stop "glr: alpha,beta > -1"

    call glr_phase_march_nodes(npts, alpha, beta, x, ders)

    wts = 1.0_dp / ((1.0_dp - x**2) * ders**2)
    C = 2**(alpha + beta + 1) * exp(log_gamma(npts + alpha + 1) - &
                                    log_gamma(npts + alpha + beta + 1) + &
                                    log_gamma(npts + beta + 1) - log_gamma(npts + 1._dp))
    wts = wts * C
end subroutine gauss_jacobi_glr

!> CAF: partition root indices; each image uses independent k-th asymptotic
!> starter + Newton polish (parallel). Under -fcoarray=single, every index is
!> owned by image 1 (same coverage as serial independent path).
subroutine gauss_jacobi_glr_caf(npts, alpha, beta, x, wts)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    real(dp), allocatable :: x_c(:)[:], d_c(:)[:]
    real(dp) :: xk, dk, C
    integer :: me, nimg, k, img

    if (npts <= 0) error stop "glr_caf: npts must be positive"
    if (alpha <= -1.0_dp .or. beta <= -1.0_dp) error stop "glr_caf: alpha,beta > -1"

    me = this_image()
    nimg = num_images()
    allocate (x_c(npts)[*], d_c(npts)[*])
    x_c = 0.0_dp
    d_c = 0.0_dp

    do k = 1, npts
        if (.not. caf_owns_index(k, me, nimg)) cycle
        call kth_root_starter(npts, alpha, beta, k, xk)
        call newton_polish(npts, alpha, beta, xk, dk)
        if (k == 1) call ensure_leftmost_root(npts, alpha, beta, xk, dk)
        x_c(k) = xk
        d_c(k) = dk
    end do
    sync all

    if (me == 1) then
        do k = 1, npts
            img = caf_root_owner(k, nimg)
            if (img == 1) then
                x(k) = x_c(k)
                xk = d_c(k)
            else
                x(k) = x_c(k)[img]
                xk = d_c(k)[img]
            end if
            ! store ders temporarily in wts, then convert
            wts(k) = xk
        end do
        x_c = x
        d_c = wts
    end if
    sync all
    x = x_c(:)[1]
    wts = d_c(:)[1]  ! ders
    sync all

    ! Weights from ders (same formula as serial glr)
    do k = 1, npts
        dk = wts(k)
        wts(k) = 1.0_dp / ((1.0_dp - x(k) * x(k)) * dk * dk)
    end do
    C = 2**(alpha + beta + 1) * exp(log_gamma(npts + alpha + 1) - &
                                    log_gamma(npts + alpha + beta + 1) + &
                                    log_gamma(npts + beta + 1) - log_gamma(npts + 1._dp))
    wts = wts * C
    deallocate (x_c, d_c)
end subroutine gauss_jacobi_glr_caf

!> Sequential Prüfer phase-by-π march: first root from endpoint asymptotic +
!> Newton; each later root from integrating ρ through Δθ = π from the previous
!> polished zero, then Newton.
subroutine glr_phase_march_nodes(npts, alpha, beta, x, ders)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), ders(npts)
    integer :: k
    real(dp) :: xk, dk, xprev, x_lo, gap

    ! Leftmost root: Bessel-endpoint asymptotic (not bulk cosθ) then Newton.
    ! Cosθ starters drift right for large α,β and Newton can latch onto root 2,
    ! which skips a zero and collapses the chain (duplicate right end).
    call first_root_starter(npts, alpha, beta, xk)
    call newton_polish(npts, alpha, beta, xk, dk)
    ! If Newton jumped right of the true first zero, re-bracket from the left.
    call ensure_leftmost_root(npts, alpha, beta, xk, dk)
    x(1) = xk
    ders(1) = dk

    do k = 2, npts
        xprev = x(k - 1)
        ! March phase by exactly π past the previous zero → next-root guess
        xk = phase_march_pi(npts, alpha, beta, xprev)
        call newton_polish(npts, alpha, beta, xk, dk)
        ! Enforce strict increase: bracket (xprev, xk] if Newton retreated
        if (xk <= xprev + 10.0_dp * x_eps) then
            gap = max(1.0e-8_dp, 2.0_dp / real(npts, dp))
            x_lo = xprev + gap * 0.25_dp
            xk = min(1.0_dp - 10.0_dp * x_eps, xprev + gap)
            call newton_bracket(npts, alpha, beta, x_lo, xk, dk)
        end if
        ! Guard against latching onto the previous root (same basin)
        if (abs(xk - xprev) < 1.0e-12_dp * (1.0_dp + abs(xprev))) then
            xk = min(1.0_dp - 10.0_dp * x_eps, xprev + max(1.0e-6_dp, pi / &
                 (rho_jacobi(npts, alpha, beta, xprev) * real(npts, dp))))
            call newton_polish(npts, alpha, beta, xk, dk)
        end if
        x(k) = xk
        ders(k) = dk
    end do
end subroutine glr_phase_march_nodes

!> First positive zero of J_ν (McMahon / NIST leading terms). ν ≥ 0.
pure real(dp) function bessel_j_first_zero(nu) result(j)
    real(dp), intent(in) :: nu
    real(dp) :: cbrt, inv

    if (nu <= 0.0_dp) then
        ! J_0 first zero
        j = 2.404825557695773_dp
        return
    end if
    ! j_{ν,1} ≈ ν + 1.8557571 ν^{1/3} + 1.033150 / ν^{1/3}  (s=1 McMahon)
    cbrt = nu**(1.0_dp / 3.0_dp)
    inv = 1.0_dp / cbrt
    j = nu + 1.855757081467_dp * cbrt + 1.033150_dp * inv
end function bessel_j_first_zero

!> Asymptotic starter for the *leftmost* Jacobi zero (near x = −1).
subroutine first_root_starter(npts, alpha, beta, x0)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x0
    call kth_root_starter(npts, alpha, beta, 1, x0)
end subroutine first_root_starter

!> Independent starter for the k-th root (k=1 leftmost … k=n rightmost).
!> Interior: cosθ bulk; near left use Bessel-β, near right Bessel-α.
subroutine kth_root_starter(npts, alpha, beta, k, x0)
    integer, intent(in) :: npts, k
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x0
    real(dp) :: nn, ab, Nscale, theta, jnu, frac

    nn = real(npts, dp)
    ab = alpha + beta
    Nscale = nn + 0.5_dp * (ab + 1.0_dp)
    frac = real(k, dp) / max(1.0_dp, nn)

    if (frac < 0.15_dp .or. (k == 1 .and. max(abs(alpha), abs(beta)) > 1.0_dp)) then
        ! Near x=-1: j_{β,s} scale for the s-th root from the left (s≈k)
        jnu = bessel_j_first_zero(max(0.0_dp, beta))
        ! crude multi-root: use McMahon shift (k-1)*pi for higher left roots
        if (k > 1) jnu = jnu + real(k - 1, dp) * pi
        x0 = -1.0_dp + (jnu * jnu) / (2.0_dp * Nscale * Nscale)
    else if (frac > 0.85_dp) then
        jnu = bessel_j_first_zero(max(0.0_dp, alpha))
        if (k < npts) jnu = jnu + real(npts - k, dp) * pi
        x0 = 1.0_dp - (jnu * jnu) / (2.0_dp * Nscale * Nscale)
    else
        ! k from left: largest θ at k=1
        theta = pi * (nn - real(k, dp) + 1.0_dp - 0.25_dp + 0.5_dp * beta) / &
                (nn + 0.5_dp * (ab + 1.0_dp))
        x0 = cos(theta)
    end if

    if (x0 <= -1.0_dp) x0 = -1.0_dp + 100.0_dp * x_eps
    if (x0 >= 1.0_dp) x0 = 1.0_dp - 100.0_dp * x_eps
end subroutine kth_root_starter

!> If the polished root is not the leftmost (P changes sign further left, or |P|
!> large while a left interval still has a zero), pull left with bisection.
subroutine ensure_leftmost_root(npts, alpha, beta, xk, dk)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(inout) :: xk
    real(dp), intent(out) :: dk
    real(dp) :: p, pp, p_left, x_left, x_try
    integer :: it

    call eval_jacobi_scalar(xk, npts, alpha, beta, p, pp)
    x_left = -1.0_dp + 100.0_dp * x_eps
    call eval_jacobi_scalar(x_left, npts, alpha, beta, p_left, pp)

    ! Sign change on (-1, xk) ⇒ a root was skipped to the left
    if (p_left * p < 0.0_dp .or. abs(p) > 1.0e-6_dp * (1.0_dp + abs(pp))) then
        if (p_left * p < 0.0_dp) then
            call newton_bracket(npts, alpha, beta, x_left, xk, dk)
            return
        end if
        ! Walk left with small steps looking for a sign change
        x_try = xk
        do it = 1, 64
            x_try = x_try - max(1.0e-4_dp, 0.5_dp / real(npts, dp))
            if (x_try <= x_left) exit
            call eval_jacobi_scalar(x_try, npts, alpha, beta, p, pp)
            if (p_left * p < 0.0_dp) then
                xk = x_try
                call newton_bracket(npts, alpha, beta, x_left, xk, dk)
                return
            end if
        end do
    end if
    call eval_jacobi_scalar(xk, npts, alpha, beta, p, pp)
    dk = pp
end subroutine ensure_leftmost_root

!> Bisection then Newton on a bracket (lo, hi] known to contain one zero.
subroutine newton_bracket(npts, alpha, beta, lo, hi, dk)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(inout) :: lo, hi
    real(dp), intent(out) :: dk
    real(dp) :: a, b, mid, pa, pb, pm, pp
    integer :: it

    a = lo
    b = hi
    if (a > b) then
        mid = a
        a = b
        b = mid
    end if
    a = clamp_x(a)
    b = clamp_x(b)
    call eval_jacobi_scalar(a, npts, alpha, beta, pa, pp)
    call eval_jacobi_scalar(b, npts, alpha, beta, pb, pp)
    if (pa * pb > 0.0_dp) then
        ! No proven bracket: Newton from midpoint
        mid = 0.5_dp * (a + b)
        call newton_polish(npts, alpha, beta, mid, dk)
        hi = mid
        return
    end if
    do it = 1, 60
        mid = 0.5_dp * (a + b)
        call eval_jacobi_scalar(mid, npts, alpha, beta, pm, pp)
        if (abs(b - a) < 1.0e-14_dp * (1.0_dp + abs(mid))) exit
        if (pa * pm <= 0.0_dp) then
            b = mid
            pb = pm
        else
            a = mid
            pa = pm
        end if
    end do
    mid = 0.5_dp * (a + b)
    call newton_polish(npts, alpha, beta, mid, dk)
    hi = mid
end subroutine newton_bracket

!> Integrate dθ/dx = ρ(x) from just right of x_prev until Δθ = π (RK4).
!> Equivalently integrate dx/dθ = 1/ρ with θ: 0 → π. Returns x_next guess.
function phase_march_pi(npts, alpha, beta, x_prev) result(x_next)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta, x_prev
    real(dp) :: x_next
    real(dp) :: th, h, xx, k1, k2, k3, k4, rho
    integer :: s, nstep

    ! Start slightly to the right of the previous zero (increasing x)
    xx = x_prev + max(1.0e-12_dp, 1.0e-4_dp / real(npts, dp))
    if (xx >= 1.0_dp - x_eps) then
        x_next = 1.0_dp - 10.0_dp * x_eps
        return
    end if

    ! More RK steps when α,β large (endpoint corrections steepen ρ)
    nstep = rk_steps_per_pi
    if (max(abs(alpha), abs(beta)) > 2.0_dp) nstep = 2 * rk_steps_per_pi
    if (max(abs(alpha), abs(beta)) > 5.0_dp) nstep = 4 * rk_steps_per_pi

    h = pi / real(nstep, dp)
    th = 0.0_dp
    do s = 1, nstep
        rho = rho_jacobi(npts, alpha, beta, xx)
        if (rho <= 0.0_dp) exit
        ! dx/dθ = 1/ρ
        k1 = 1.0_dp / rho
        rho = rho_jacobi(npts, alpha, beta, clamp_x(xx + 0.5_dp * h * k1))
        k2 = 1.0_dp / max(rho, epsilon(1.0_dp))
        rho = rho_jacobi(npts, alpha, beta, clamp_x(xx + 0.5_dp * h * k2))
        k3 = 1.0_dp / max(rho, epsilon(1.0_dp))
        rho = rho_jacobi(npts, alpha, beta, clamp_x(xx + h * k3))
        k4 = 1.0_dp / max(rho, epsilon(1.0_dp))
        xx = clamp_x(xx + (h / 6.0_dp) * (k1 + 2.0_dp * k2 + 2.0_dp * k3 + k4))
        th = th + h
        if (xx >= 1.0_dp - x_eps) exit
    end do
    x_next = xx
end function phase_march_pi

!> Local wave number for Jacobi Prüfer (bulk LG + endpoint corrections).
!> ρ² ≈ λ/(1−x²) + (1−α²)/(4(1−x)²) + (1−β²)/(4(1+x)²), λ = n(n+α+β+1).
!> When endpoint terms make r2 negative (α or β > 1 near the ends), floor with
!> a fraction of the bulk term so phase march still advances.
pure real(dp) function rho_jacobi(n, alpha, beta, x) result(rho)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha, beta, x
    real(dp) :: lam, den, xm, xp, r2, bulk

    lam = real(n, dp) * (real(n, dp) + alpha + beta + 1.0_dp)
    den = 1.0_dp - x * x
    if (den < 1.0e-30_dp) den = 1.0e-30_dp
    xm = max(1.0e-14_dp, 1.0_dp - x)
    xp = max(1.0e-14_dp, 1.0_dp + x)
    bulk = lam / den
    r2 = bulk + (1.0_dp - alpha * alpha) / (4.0_dp * xm * xm) + &
         (1.0_dp - beta * beta) / (4.0_dp * xp * xp)
    if (r2 < 0.05_dp * bulk) r2 = 0.05_dp * bulk
    if (r2 < 0.0_dp) r2 = bulk
    rho = sqrt(r2)
end function rho_jacobi

pure real(dp) function clamp_x(x) result(y)
    real(dp), intent(in) :: x
    y = x
    if (y <= -1.0_dp) y = -1.0_dp + 10.0_dp * x_eps
    if (y >= 1.0_dp) y = 1.0_dp - 10.0_dp * x_eps
end function clamp_x

subroutine newton_polish(npts, alpha, beta, xk, dk)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(inout) :: xk
    real(dp), intent(out) :: dk
    real(dp) :: p, pp, dx, x_old
    integer :: it

    do it = 1, newton_max_it
        call eval_jacobi_scalar(xk, npts, alpha, beta, p, pp)
        if (abs(pp) < epsilon(1.0_dp)) exit
        dx = -p / pp
        ! Damp large steps so we do not jump over a neighboring zero
        if (abs(dx) > 0.25_dp) dx = sign(0.25_dp, dx)
        x_old = xk
        xk = clamp_x(xk + dx)
        if (abs(dx) < 1.0e-14_dp * (1.0_dp + abs(xk))) exit
        if (abs(xk - x_old) < 1.0e-16_dp) exit
    end do
    call eval_jacobi_scalar(xk, npts, alpha, beta, p, pp)
    dk = pp
end subroutine newton_polish

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

end module gjp_glr

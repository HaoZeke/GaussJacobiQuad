!> @brief This module provides routines for numerical integration using Gauss-Jacobi quadrature.
!! @details For degree less than 4, the implementation uses the standard Golub-Welsch method.
!! For 4 or more points, the implementation adopts the Gil-Segura-Temme (GST)
!! method as described in their 2021 paper and accompanying (MIT licensed) MATLAB code.
!!
!! References:
!!  - Gil, A., Segura, J. & Temme, N.M. (2021). Fast and reliable high-accuracy computation of Gauss–Jacobi quadrature.
!!    Numerical Algorithms, 87, 1391–1419. DOI: 10.1007/s11075-020-01012-6
!!  - Golub, G. H., & Welsch, J. H. (1969). Calculation of Gauss Quadrature Rules. Mathematics of Computation, 23(106), 221–230.
!!    DOI: 10.2307/2004418
!!
module gjp_gst
use gjp_types, only: dp
use gjp_common, only: jacobi_matrix, jacobi_zeroeth_moment
implicit none
real(dp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510
contains

subroutine gauss_jacobi_gst(npts, alpha, beta, x, wts)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    ! Local variables
    real(dp) :: lo, s
    integer :: n_x1, n_x2, idx, total_length
    real(dp), allocatable :: x1(:), x2(:), w1(:), w2(:), w0(:), x0(:)
    integer :: jl, jr, numec, le, i, n
    real(dp) :: t, tn, wn, del, del1, del2, la0, la1, la2
    real(dp) :: la, r1, r2, s00, s01, s02, s10, s11, s12, s20, s21, s22
    real(dp) :: a, b
    a = alpha
    b = beta
    n = npts

    lo = (alpha + beta + 1.0_dp) * log(2.0_dp) &
         + log_gamma(alpha + 1.0_dp) &
         + log_gamma(beta + 1.0_dp) &
         - log_gamma(alpha + beta + 2.0_dp)
    ! fsweep populates x1 and w1 arrays.
    call fsweep(npts, alpha, beta, x1, w1)

    ! fsweep populates x2 and w2 arrays.
    call fsweep(npts, beta, alpha, x2, w2)
    s = 0.0_dp

    ! Check if x1 is empty
    if (allocated(x1)) then
        n_x1 = size(x1)
        if (x1(1) < 0.0_dp) then
            ! Remove first element
            x1 = x1(2:n_x1)
            w1 = w1(2:n_x1)
        end if
    end if
    if (allocated(x2)) then
        n_x2 = size(x2)
        if (x2(1) < 0.0_dp) then
            x2 = x2(2:n_x2)
            w2 = w2(2:n_x2)
        end if
    end if
    if (allocated(x1) .and. allocated(x2)) then
        ! ?? Magic number?
        if (abs(x2(1) + x(1)) < 3.0e-16_dp) then
            s = 1.0_dp
        end if
    end if

    ! Reverse x2 and w2 and store in x and w
    do idx = 1, n_x2 - s
        x(idx) = x2(n_x2 + 1 - idx - s)
        wts(idx) = w2(n_x2 + 1 - idx - s)
    end do
    ! Append x1 and w1 to x and w
    do idx = 1, n_x1
        x(n_x2 - s + idx) = x1(idx)
        wts(n_x2 - s + idx) = w1(idx)
    end do

    ! Compute fixed point
    jl = 0
    jr = 0
    ! Calculate numec
    numec = max(3, int(log(real(npts, dp)) / log(10.0_dp)) + 1)

    do idx = 1, numec
        if ((x(idx) < -0.95_dp) .and. (beta < -0.01_dp)) then
            jl = jl + 1
            t = acos(-x(idx))
            call fixedt(npts, beta, alpha, tn, wn)
            wts(idx) = wn
            x(idx) = -cos(tn)
        end if
    end do
    le = size(x)
    do idx = le, le - numec + 1, -1
        if ((x(idx) > 0.95_dp) .and. (alpha < -0.01_dp)) then
            jr = jr + 1
            t = acos(x(idx))
            ! Here, fixedt needs to be a subroutine that modifies tn and wn
            call fixedt(npts, alpha, beta, t, wn)
            wts(idx) = wn
            x(idx) = cos(tn)
        end if
    end do
    ! Partition the w array into w0, w1, w2
    w0 = wts(1:jl)
    w1 = wts(jl + 1:le - jr)
    w2 = wts(le - jr + 1:le)
    ! Normalize w1
    w1 = w1 / w1(1)
    total_length = size(w0) + size(w1) + size(w2)
    ! Calculate w based on conditions
    if (a > -0.75_dp .and. b > -0.75_dp) then
        w0 = fa(n, b, a) * w0
        w2 = fa(n, a, b) * w2
        la = (exp(lo) - sum(w0) - sum(w2)) / sum(w1)
        wts(1:size(w0)) = w0
        wts(size(w0) + 1:size(w0) + size(w1)) = la * w1
        wts(size(w0) + size(w1) + 1:total_length) = w2
    else
        ! Partition the x array into x0, x1, x2

        allocate (x0(jl))
        allocate (x1(le - jr - jl + 1))
        allocate (x2(jr))
        do i = 1, jl
            x0(i) = x(i)
        end do
        do i = 1, le - jr - jl + 1
            x1(i) = x(jl + i)
        end do

        do i = 1, jr
            x2(i) = x(le - jr + i)
        end do

        ! Compute other variables
        r1 = (b - a) / (a + b + 2)
        r2 = ((a - b)**2 + a + b + 2) / ((a + b + 2) * (a + b + 3))

        ! Compute sums
        s00 = sum(w0)
        s01 = sum(w1)
        s02 = sum(w2)

        s10 = sum(x0 * w0)
        s11 = sum(x1 * w1)
        s12 = sum(x2 * w2)

        s20 = sum(x0**2 * w0)
        s21 = sum(x1**2 * w1)
        s22 = sum(x2**2 * w2)

        if (s00 * s01 * s02 > 0) then
            del1 = s00 * s11 * s22 + s20 * s01 * s12 + s10 * s21 * s02
            del2 = s20 * s11 * s02 + s10 * s01 * s22 + s00 * s12 * s21
            del = del1 - del2
            la0 = s11 * s22 + r2 * s01 * s12 + r1 * s02 * s21 - s21 * s12 - r1 * s01 * s22 - r2 * s11 * s02
            la0 = la0 / del
            la1 = -s10 * s22 - r2 * s00 * s12 - r1 * s02 * s20 + s20 * s12 + r1 * s00 * s22 + r2 * s10 * s02
            la1 = la1 / del
            la2 = s10 * s21 + r2 * s00 * s11 + r1 * s01 * s20 - s20 * s11 - r1 * s00 * s21 - r2 * s10 * s01
            la2 = la2 / del
            wts(1:jl) = la0 * w0
            wts(jl + 1:le - jr) = la1 * w1
            wts(le - jr + 1:le) = la2 * w2
        else if (s00 * s01 > 0) then
            del = s00 * s11 - s10 * s01
            la0 = (s11 - r1 * s01) / del
            la1 = (r1 * s00 - s10) / del
            wts(1:jl) = la0 * w0
            wts(jl:size(wts)) = la1 * w1
        else if (s01 * s02 > 0) then
            del = s01 * s12 - s11 * s02
            la1 = (s12 - r1 * s02) / del
            la2 = (r1 * s01 - s11) / del
            wts(1:jl) = la0 * w0
            wts(jl:size(wts)) = la2 * w2
        else
            wts = wts / sum(wts)
        end if
        wts = wts * exp(lo)
    end if
end subroutine gauss_jacobi_gst

subroutine fsweep(n, a, b, x, w)
    ! Amatrgument declarations
    integer, intent(in) :: n
    real(dp), intent(in) :: a, b
    real(dp), allocatable, intent(out) :: x(:), w(:)

    ! Variable declarations
    real(dp) :: L, Bmat, Amat, Delta, xf, xi, nor, facm
    real(dp) :: xt, wt, y0, y1, sqde, de, t, alpha, beta
    integer :: i, index
    ! Initialize variables
    nor = 1.0d-130
    L = 2 * n + a + b + 1
    Bmat = (b - a) * ((a + b + 6) * n + 2 * (a + b))
    Amat = (2 * n + a + b) * (n * (2 * n + a + b) + 2 * (a + b + 2))
    Delta = n**2 * (n + a + b + 1)**2 + (a + 1) * (b + 1) * &
            (n**2 + (a + b + 4) * n + 2 * (a + b))
    xf = (Bmat + 4 * (n - 1) * sqrt(Delta)) / Amat
    xi = (b**2 - a**2) / (L**2 - 1)
    alpha = a
    beta = b

    ! Initialize y0, y1
    if ((beta - alpha) == 0.0_dp) then
        index = n - 2 * floor(real(n / 2), dp)
        if (index > 0) then
            y0 = 0.0_dp
            y1 = nor
        else
            y0 = nor
            y1 = 0.0_dp
        end if
    else
        y1 = nor
        facm = recu(n, alpha, beta, xi)
        if (abs(facm) > 100) then
            xi = 2 * xi
            facm = recu(n, alpha, beta, xi)
        end if
        y0 = y1 * facm
    end if
    ! Initialize array sizes
    allocate (x(n))
    allocate (w(n))

    ! Main loop
    i = 1
    do while (xi < xf)
        ! Calling fixedzser subroutine (insert your actual implementation)
        call fixedzser(n, a, b, xi, xt, wt, y0, y1)

        ! Increment index and set values
        i = i + 1
        x(i) = xt
        w(i) = wt

        ! Check for termination condition
        if (i > 3 .and. abs(x(i) / x(i - 1) - 1) < 2.d-15) then
            xf = 0
            x = x(1:i - 2)
            w = w(1:i - 2)
            exit
        end if

        ! Compute delta
        sqde = sqrt(0.25d0 * ((L * L - 1) * (1 - xt * xt) - 2 * a * a * (1 + xt) - 2 * b * b * (1 - xt)))
        de = -tanh(pi / sqde)
        xi = (xt - de) / (1 - xt * de)
        t = xi - xt

        ! Update y0, y1
        ! Calling seriestay subroutine (insert your actual implementation)
        call seriestay(n, a, b, xt, t, y0, y1)

    end do

    ! Resize arrays to actual size
    if (i < n) then
        x = x(1:i)
        w = w(1:i)
    end if
end subroutine fsweep

function recu(n, a, b, x) result(r)
    integer, intent(in) :: n
    real(dp), intent(in) :: a, b, x
    real(dp) :: r ! Y/Y' for initial Taylor expansion
    real(dp) :: h, L, Amat, Bmat, C
    ! h is P(n, a, b, x) / P(n - 1, a, b, x)
    integer :: k

    h = 0.5_dp * (a - b + (a + b + 2.0_dp) * x)

    do k = 1, n - 1
        L = 2.0_dp * k + a + b + 1.0_dp
        Amat = 2.0_dp * (k + 1) * (k + a + b + 1) * (L - 1.0_dp)
        Bmat = L * ((L**2 - 1.0_dp) * x + a**2 - b**2)
        C = 2.0_dp * (k + a) * (k + b) * (L + 1.0_dp)
        h = (Bmat - C / h) / Amat
    end do

    L = 2.0_dp * n + a + b + 1.0_dp
    r = (n + b + 1.0_dp) / 2.0_dp / (1.0_dp + x) &
        - (n + a + 1.0_dp) / 2.0_dp / (1.0_dp - x) &
        + 1.0_dp / (L - 1.0_dp) / (1.0_dp - x**2) &
        * (n * (a - b) + 2.0_dp * (n + a) * (n + b) / h)
    r = 1.0_dp / r

end function recu

subroutine fixedzser(n, a, b, x0, y0, y1, x, w)
! fixed point on z and using Taylor for (1-x)^(a+1)/2*(1+x)^(b+1)/2*P
    integer, intent(in) :: n
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: x, w
    real(dp), intent(inout) :: x0, y0, y1
    real(dp) :: epsil, erro, L, sq, sqde, h, argu, de, t
    integer :: i

    epsil = 2.3e-16_dp
    erro = 1.0_dp
    L = 2.0_dp * n + a + b + 1.0_dp
    i = 0
    x = x0

    do while (erro > epsil)
        i = i + 1
        sq = 0.25_dp * ((L * L - 1.0_dp) * (1.0_dp - x * x) - 2.0_dp * a * a * (1.0_dp + x) - 2.0_dp * b * b * (1.0_dp - x))
        sqde = sqrt(abs(sq))
        h = y0 / (y1 * (1.0_dp - x * x) + x * y0)

        if (sq > 0.0_dp) then
            argu = atan(sqde * h) / sqde
            if (h > 1.0_dp) then
                argu = argu - pi / sqde
            end if
            de = tanh(argu)
        else
            de = tanh(atanh(sqde * h) / sqde)
        end if

        x = (x - de) / (1.0_dp - x * de)
        erro = abs(1.0_dp - x / x0)
        t = x - x0
        call seriestay(n, a, b, x0, t, y0, y1)
        x0 = x
    end do

    w = (1.0_dp - x)**a * (1.0_dp + x)**b / y1 / y1

end subroutine fixedzser

subroutine seriestay(n, a, b, xp, t, y0, yp1)
    integer, intent(in) :: n
    real(dp), intent(in) :: a, b, xp, t
    real(dp), intent(inout) :: y0, yp1
    real(dp) :: epso, delta, suma, sumad, errod, xp2, a2, b2, L, L2
    real(dp) :: ym1, ym2, yp2, nt, ntd
    integer :: j
    real(dp) :: c1, c2, c3, c4, c5, tx

    epso = 1.0e-19_dp
    delta = t
    ym1 = 0.0_dp
    ym2 = 0.0_dp
    y0 = y0 * 1.0e-30_dp
    yp1 = yp1 * 1.0e-30_dp
    suma = y0 + yp1 * t
    sumad = yp1
    j = -1
    errod = 1.0_dp
    xp2 = xp * xp
    a2 = a * a
    b2 = b * b
    L = 2.0_dp * n + a + b + 1.0_dp
    L2 = L * L

    do while (errod > epso .and. j < 100)
        j = j + 1
        c1 = 4.0_dp * (j + 2.0_dp) * (j + 1) * (1.0_dp - xp2)**2
        c2 = -16.0_dp * (j + 1) * j * (1 - xp2) * xp
        c3 = 0.5_dp * j * (j - 1) * (48 * xp2 - 16) + (L2 - 1) * (1 - xp2) &
             - 2 * (b2 - 1) * (1 - xp) - 2 * (a2 - 1) * (1 + xp)
        c4 = 16 * (j - 1) * (j - 2) * xp - 2 * (L2 - 1) * xp + 2 * b2 - 2 * a2
        c5 = 4 * (j - 2) * (j - 3) + (1 - L2)
        yp2 = -1.0_dp / c1 * (c2 * yp1 + c3 * y0 + c4 * ym1 + c5 * ym2)

        ym2 = ym1
        ym1 = y0
        y0 = yp1
        yp1 = yp2

        ntd = (j + 2) * yp2 * t
        sumad = sumad + ntd
        tx = t * delta
        nt = yp2 * tx
        suma = suma + nt

        if (j > 20.0_dp .and. y0 * yp1 /= 0.0_dp) then
            errod = max(abs(ntd / sumad), abs(nt / suma))
        end if
    end do

    y0 = suma * 1.0e30_dp
    yp1 = sumad * 1.0e30_dp

end subroutine seriestay

subroutine fixedt(n, a, b, t, w)
    integer, intent(in) :: n
    real(dp), intent(in) :: a, b
    real(dp), intent(inout) :: t
    real(dp), intent(out) :: w
    real(dp) :: epsil, erro, L, st2, sq, sqde, h, at, t0
    integer :: i

    epsil = 6.0e-16_dp
    erro = 1.0_dp
    L = 2.0_dp * n + a + b + 1.0_dp
    i = 0
    t0 = t

    do while (erro > epsil)
        i = i + 1
        st2 = sin(t0 / 2.0_dp)**2
        sq = 0.25_dp - a**2 + 0.25_dp * L**2 * sin(t0)**2 + (a**2 - b**2) * st2
        sqde = sqrt(abs(sq))
        h = 0.5_dp + a + L * st2 - 2.0_dp * (n + a + b + 1) * st2 * cft(n, a + 1.0_dp, b, t0)
        h = 1.0_dp / h

        if (sq > 0.0_dp) then
            at = atan(sqde * h)
        else
            at = atanh(sqde * h)
        end if

        t = t0 - sin(t0) / sqde * at
        erro = abs(1.0_dp - t / t0)
        t0 = t
    end do

    ! Assuming seriest is another subroutine you've implemented that takes the same parameters.
    ! Replace with the actual function if different.
    w = 1.0_dp / seriest(n - 1, a + 1.0_dp, b + 1.0_dp, t)**2 / sin(t)**2
end subroutine fixedt

! Implementation of cft function
function cft(n, a, b, theta) result(cf)
    integer, intent(in) :: n
    real(dp), intent(in) :: a, b, theta
    real(dp) :: cf
    real(dp) :: ai, t, bp, ap, bb, aa, E, F, del, R

    ai = a + 1.0_dp
    t = sin(theta / 2.0_dp)**2
    bp = -1.0_dp - (n * t + a) / ((n + a + b + 1) * t)
    ap = -(a + n) / (n + a + b + 1) / t
    bb = -1.0_dp - (n * t + ai) / ((n + ai + b + 1) * t)
    aa = -(ai + n) / (n + ai + b + 1) / t

    E = bb
    F = bb + aa / bp
    cf = ap / (bp + aa / bb)
    del = abs(E / F - 1.0_dp)

    do while (del > 1.0e-15_dp)
        ai = ai + 1.0_dp
        bb = -1.0_dp - (n * t + ai) / ((n + ai + b + 1) * t)
        aa = -(ai + n) / (n + ai + b + 1) / t
        E = bb + aa / E
        F = bb + aa / F
        R = E / F
        del = abs(R - 1.0_dp)
        cf = cf * R
    end do
end function cft

! Implementation of seriest function
function seriest(n, a, b, theta) result(s)
    integer, intent(in) :: n
    real(dp), intent(in) :: a, b, theta
    real(dp) :: s
    real(dp) :: t, de, er
    integer :: k

    t = sin(theta / 2.0_dp)**2
    s = 1.0_dp
    de = 1.0_dp
    k = 0
    er = 1.0_dp

    do while (er > 1.0e-17_dp)
        k = k + 1
        de = (k + n + a + b) * (k - 1 - n) / (k * (k + a)) * t * de
        s = s + de
        er = abs(de / s)
    end do
end function seriest

function fa(n, a, b) result(fa_value)
    ! Argument declarations
    integer, intent(in) :: n
    real(dp), intent(in) :: a, b
    real(dp) :: fa_value ! The result

    ! Local variables
    real(dp) :: lofa1, lofa2

    ! Compute lofa1 and lofa2 based on given formula
    lofa1 = (a + b + 3.0_dp) * log(2.0_dp) + 2.0_dp * log_gamma(a + 2.0_dp) + &
            log_gamma(real(n, dp)) + log_gamma(real(n, dp) + b + 1.0_dp)
    lofa2 = log(real(n, dp)) + &
            2.0_dp * log(real(n, dp) + a + b + 1.0_dp) + &
            log_gamma(real(n, dp) + a + 1.0_dp) + &
            log_gamma(real(n, dp) + a + b + 1.0_dp)

    ! Compute fa_value
    fa_value = exp(lofa1 - lofa2)
end function fa
end module gjp_gst

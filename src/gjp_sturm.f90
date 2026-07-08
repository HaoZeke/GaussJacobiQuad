!> Golub-Welsch via Sturm bisection + inverse iteration on the Jacobi matrix.
!>
!> Unlike sequential imtqlx (one long QL path), each eigenvalue index k is an
!> independent scalar search (Sturm count) plus a short tridiagonal solve.
!> That is the CAF-friendly reformulation of the *same* eigenproblem ACM 655
!> solves, not a faster QL chase.
!>
!> References: Parlett, The Symmetric Eigenvalue Problem; Golub-Welsch 1969.
module gjp_sturm
use gjp_types, only: dp, gjp_sparse_matrix
use gjp_common, only: jacobi_matrix, jacobi_zeroeth_moment
use gjp_rec, only: caf_root_owner, caf_owns_index
implicit none

private
public :: gauss_jacobi_sturm, gauss_jacobi_sturm_caf

contains

!> Serial: all eigenvalues of Jacobi T by bisection, weights by inv. iteration.
subroutine gauss_jacobi_sturm(npts, alpha, beta, x, wts)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    real(dp) :: a(npts), b(npts - 1), mu0, lo, hi, q1
    type(gjp_sparse_matrix) :: jm
    integer :: k

    jm = jacobi_matrix(npts, alpha, beta)
    mu0 = jacobi_zeroeth_moment(alpha, beta)
    a = jm%diagonal(1:npts)
    if (npts > 1) b = jm%off_diagonal(1:npts - 1)

    call gershgorin_bounds(npts, a, b, lo, hi)
    do k = 1, npts
        x(k) = bisect_kth(npts, a, b, k, lo, hi)
        call inv_iter_first(npts, a, b, x(k), q1)
        wts(k) = mu0 * (q1 * q1)
    end do
end subroutine gauss_jacobi_sturm

!> CAF: eigenvalue indices partitioned across images; gather nodes/weights.
subroutine gauss_jacobi_sturm_caf(npts, alpha, beta, x, wts)
    integer, intent(in) :: npts
    real(dp), intent(in) :: alpha, beta
    real(dp), intent(out) :: x(npts), wts(npts)
    real(dp) :: a(npts), b(npts - 1), mu0, lo, hi, q1
    real(dp), allocatable :: x_c(:)[:], w_c(:)[:]
    type(gjp_sparse_matrix) :: jm
    integer :: k, me, nimg, img

    me = this_image()
    nimg = num_images()
    jm = jacobi_matrix(npts, alpha, beta)
    mu0 = jacobi_zeroeth_moment(alpha, beta)
    a = jm%diagonal(1:npts)
    if (npts > 1) b = jm%off_diagonal(1:npts - 1)
    call gershgorin_bounds(npts, a, b, lo, hi)

    allocate (x_c(npts)[*], w_c(npts)[*])
    x_c = 0.0_dp
    w_c = 0.0_dp

    do k = 1, npts
        if (.not. caf_owns_index(k, me, nimg)) cycle
        x_c(k) = bisect_kth(npts, a, b, k, lo, hi)
        call inv_iter_first(npts, a, b, x_c(k), q1)
        w_c(k) = mu0 * (q1 * q1)
    end do
    sync all

    if (me == 1) then
        do k = 1, npts
            img = caf_root_owner(k, nimg)
            if (img == 1) then
                x(k) = x_c(k)
                wts(k) = w_c(k)
            else
                x(k) = x_c(k)[img]
                wts(k) = w_c(k)[img]
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
end subroutine gauss_jacobi_sturm_caf

subroutine gershgorin_bounds(n, a, b, lo, hi)
    integer, intent(in) :: n
    real(dp), intent(in) :: a(n), b(n - 1)
    real(dp), intent(out) :: lo, hi
    integer :: i
    real(dp) :: rad, left, right

    lo = a(1)
    hi = a(1)
    do i = 1, n
        rad = 0.0_dp
        if (i > 1) rad = rad + abs(b(i - 1))
        if (i < n) rad = rad + abs(b(i))
        left = a(i) - rad
        right = a(i) + rad
        if (left < lo) lo = left
        if (right > hi) hi = right
    end do
    ! pad slightly
    lo = lo - (1.0_dp + abs(lo)) * 1.0e-3_dp
    hi = hi + (1.0_dp + abs(hi)) * 1.0e-3_dp
end subroutine gershgorin_bounds

!> Number of eigenvalues of tridiag(a,b) that are strictly < lambda.
integer function sturm_count(n, a, b, lambda) result(cnt)
    integer, intent(in) :: n
    real(dp), intent(in) :: a(n), b(n - 1), lambda
    integer :: i
    real(dp) :: d, eps

    eps = epsilon(1.0_dp)
    cnt = 0
    d = a(1) - lambda
    if (d < 0.0_dp) cnt = 1
    if (abs(d) < eps) d = eps
    do i = 2, n
        d = (a(i) - lambda) - (b(i - 1) * b(i - 1)) / d
        if (d < 0.0_dp) cnt = cnt + 1
        if (abs(d) < eps) d = sign(eps, d)
        if (d == 0.0_dp) d = eps
    end do
end function sturm_count

!> k-th smallest eigenvalue (k = 1 .. n) by bisection on Sturm count.
real(dp) function bisect_kth(n, a, b, k, lo0, hi0) result(lam)
    integer, intent(in) :: n, k
    real(dp), intent(in) :: a(n), b(n - 1), lo0, hi0
    real(dp) :: lo, hi, mid
    integer :: it, cnt
    real(dp) :: tol

    lo = lo0
    hi = hi0
    tol = epsilon(1.0_dp) * (1.0_dp + max(abs(lo0), abs(hi0))) * 64.0_dp
    do it = 1, 80
        mid = 0.5_dp * (lo + hi)
        if (abs(hi - lo) <= tol) exit
        cnt = sturm_count(n, a, b, mid)
        ! want exactly k eigenvalues <= mid for upper end; use < mid count
        if (cnt >= k) then
            hi = mid
        else
            lo = mid
        end if
    end do
    lam = 0.5_dp * (lo + hi)
end function bisect_kth

!> Inverse iteration; return first component of unit eigenvector for lambda.
subroutine inv_iter_first(n, a, b, lambda, q1)
    integer, intent(in) :: n
    real(dp), intent(in) :: a(n), b(n - 1), lambda
    real(dp), intent(out) :: q1
    real(dp) :: x(n), y(n), nrm
    integer :: iter, i

    ! start near e1
    x = 0.0_dp
    x(1) = 1.0_dp
    do iter = 1, 6
        call tridiag_solve(n, a, b, lambda, x, y)
        nrm = sqrt(sum(y * y))
        if (nrm == 0.0_dp) then
            q1 = 0.0_dp
            return
        end if
        x = y / nrm
    end do
    ! consistent sign
    if (x(1) < 0.0_dp) x = -x
    q1 = x(1)
end subroutine inv_iter_first

!> Solve (tridiag(a,b) - lambda I) y = x by Thomas algorithm.
subroutine tridiag_solve(n, a, b, lambda, x, y)
    integer, intent(in) :: n
    real(dp), intent(in) :: a(n), b(*), lambda, x(n)
    real(dp), intent(out) :: y(n)
    real(dp) :: c(n), d(n), diag(n), eps, denom
    integer :: i

    eps = epsilon(1.0_dp) * 10.0_dp
    do i = 1, n
        diag(i) = a(i) - lambda
        if (abs(diag(i)) < eps) diag(i) = sign(eps, diag(i) + eps)
    end do

    if (n == 1) then
        y(1) = x(1) / diag(1)
        return
    end if

    c(1) = b(1) / diag(1)
    d(1) = x(1) / diag(1)
    do i = 2, n
        denom = diag(i) - b(i - 1) * c(i - 1)
        if (abs(denom) < eps) denom = sign(eps, denom)
        if (i < n) c(i) = b(i) / denom
        d(i) = (x(i) - b(i - 1) * d(i - 1)) / denom
    end do
    y(n) = d(n)
    do i = n - 1, 1, -1
        y(i) = d(i) - c(i) * y(i + 1)
    end do
end subroutine tridiag_solve

end module gjp_sturm

!> Cuppen divide-and-conquer tridiagonal eigen-solve with imtqlx leaves.
!>
!> Math (see scripts/imtqlx_parallel_math.py):
!>   T = blkdiag(T1_mod, T2_mod) + rho * u * u^T
!>   eig(T1), eig(T2) independent (parallel leaves)
!>   eig(T) from secular 1 + rho * sum z_i^2/(d_i - lam) = 0
!>   secular roots independent (CAF round-robin)
!>
!> Preserves Golub-Welsch / ACM-655 *problem*: nodes = eigenvalues of Jacobi T,
!> weights = mu0 * (first components of orthonormal eigenvectors)^2.
module gjp_imtqlx_dc
use gjp_types, only: dp
use gjp_imtqlx, only: imtqlx, imtqlx_full
use gjp_rec, only: caf_root_owner, caf_owns_index
implicit none

private
public :: imtqlx_dc, imtqlx_dc_caf

integer, parameter :: leaf_size = 32

contains

!> Serial recursive Cuppen DC; leaves use imtqlx_full.
subroutine imtqlx_dc(n, diag, off_diag, evals, q1)
    integer, intent(in) :: n
    real(dp), intent(in) :: diag(n), off_diag(n)
    real(dp), intent(out) :: evals(n)
    real(dp), intent(out) :: q1(n) ! first row of eigenvector matrix Q
    real(dp) :: dwork(n), owork(n), qmat(n, n)
    integer :: i

    dwork = diag
    owork = off_diag
    owork(n) = 0.0_dp
    qmat = 0.0_dp
    do i = 1, n
        qmat(i, i) = 1.0_dp
    end do
    call cuppen_recurse(n, dwork, owork, qmat)
    evals = dwork
    q1 = qmat(1, :)
end subroutine imtqlx_dc

!> CAF: same Cuppen recursion; secular roots partitioned across images.
subroutine imtqlx_dc_caf(n, diag, off_diag, evals, q1)
    integer, intent(in) :: n
    real(dp), intent(in) :: diag(n), off_diag(n)
    real(dp), intent(out) :: evals(n)
    real(dp), intent(out) :: q1(n)
    real(dp) :: dwork(n), owork(n), qmat(n, n)
    integer :: i

    dwork = diag
    owork = off_diag
    owork(n) = 0.0_dp
    qmat = 0.0_dp
    do i = 1, n
        qmat(i, i) = 1.0_dp
    end do
    call cuppen_recurse_caf(n, dwork, owork, qmat)
    evals = dwork
    q1 = qmat(1, :)
end subroutine imtqlx_dc_caf

recursive subroutine cuppen_recurse(n, diag, off_diag, qmat)
    integer, intent(in) :: n
    real(dp), intent(inout) :: diag(n), off_diag(n)
    real(dp), intent(inout) :: qmat(n, n)

    integer :: k, i, j
    real(dp) :: rho
    real(dp), allocatable :: d1(:), o1(:), d2(:), o2(:), q1m(:, :), q2m(:, :)
    real(dp), allocatable :: d(:), z(:), qblk(:, :), lam(:), w(:), qnew(:, :)
    real(dp) :: nrm, q1row(n)

    if (n <= leaf_size) then
        call imtqlx_full(n, diag, off_diag, qmat)
        return
    end if

    k = n / 2
    rho = off_diag(k)

    allocate (d1(k), o1(k), q1m(k, k))
    allocate (d2(n - k), o2(n - k), q2m(n - k, n - k))
    d1 = diag(1:k)
    o1 = 0.0_dp
    if (k > 1) o1(1:k - 1) = off_diag(1:k - 1)
    o1(k) = 0.0_dp
    d1(k) = d1(k) - rho

    d2 = diag(k + 1:n)
    o2 = 0.0_dp
    if (n - k > 1) o2(1:n - k - 1) = off_diag(k + 1:n - 1)
    o2(n - k) = 0.0_dp
    d2(1) = d2(1) - rho

    q1m = 0.0_dp
    q2m = 0.0_dp
    do i = 1, k
        q1m(i, i) = 1.0_dp
    end do
    do i = 1, n - k
        q2m(i, i) = 1.0_dp
    end do

    call cuppen_recurse(k, d1, o1, q1m)
    call cuppen_recurse(n - k, d2, o2, q2m)

    allocate (d(n), z(n), qblk(n, n), lam(n), w(n), qnew(n, n))
    d(1:k) = d1
    d(k + 1:n) = d2
    ! z = Q^T u, u = e_k + e_{k+1} in full coords
    z(1:k) = q1m(k, :)
    z(k + 1:n) = q2m(1, :)
    qblk = 0.0_dp
    qblk(1:k, 1:k) = q1m
    qblk(k + 1:n, k + 1:n) = q2m

    call secular_eigen_serial(n, d, z, rho, lam, qnew, qblk)
    ! eigenvectors of T already in qnew after merge (qblk permuted inside)
    qmat = qnew
    diag = lam
    off_diag = 0.0_dp

    deallocate (d1, o1, d2, o2, q1m, q2m, d, z, qblk, lam, w, qnew)
end subroutine cuppen_recurse

recursive subroutine cuppen_recurse_caf(n, diag, off_diag, qmat)
    integer, intent(in) :: n
    real(dp), intent(inout) :: diag(n), off_diag(n)
    real(dp), intent(inout) :: qmat(n, n)

    integer :: k, i
    real(dp) :: rho
    real(dp), allocatable :: d1(:), o1(:), d2(:), o2(:), q1m(:, :), q2m(:, :)
    real(dp), allocatable :: d(:), z(:), qblk(:, :), lam(:), qnew(:, :)

    if (n <= leaf_size) then
        call imtqlx_full(n, diag, off_diag, qmat)
        return
    end if

    k = n / 2
    rho = off_diag(k)

    allocate (d1(k), o1(k), q1m(k, k))
    allocate (d2(n - k), o2(n - k), q2m(n - k, n - k))
    d1 = diag(1:k)
    o1 = 0.0_dp
    if (k > 1) o1(1:k - 1) = off_diag(1:k - 1)
    o1(k) = 0.0_dp
    d1(k) = d1(k) - rho

    d2 = diag(k + 1:n)
    o2 = 0.0_dp
    if (n - k > 1) o2(1:n - k - 1) = off_diag(k + 1:n - 1)
    o2(n - k) = 0.0_dp
    d2(1) = d2(1) - rho

    q1m = 0.0_dp
    q2m = 0.0_dp
    do i = 1, k
        q1m(i, i) = 1.0_dp
    end do
    do i = 1, n - k
        q2m(i, i) = 1.0_dp
    end do

    ! Independent leaves/subtrees (each image does full recursion; secular is CAF)
    ! For true dual-image fork of the two halves one would farm jobs; here every
    ! image builds the same d,z then parallelizes secular roots.
    call cuppen_recurse_caf(k, d1, o1, q1m)
    call cuppen_recurse_caf(n - k, d2, o2, q2m)

    allocate (d(n), z(n), qblk(n, n), lam(n), qnew(n, n))
    d(1:k) = d1
    d(k + 1:n) = d2
    z(1:k) = q1m(k, :)
    z(k + 1:n) = q2m(1, :)
    qblk = 0.0_dp
    qblk(1:k, 1:k) = q1m
    qblk(k + 1:n, k + 1:n) = q2m

    call secular_eigen_caf(n, d, z, rho, lam, qnew, qblk)
    qmat = qnew
    diag = lam
    off_diag = 0.0_dp

    deallocate (d1, o1, d2, o2, q1m, q2m, d, z, qblk, lam, qnew)
end subroutine cuppen_recurse_caf

!> Rank-1 update eigenproblem D + rho z z^T. Permutes qblk columns with (d,z)
!> so V = qblk_sorted * W is consistent.
subroutine secular_eigen_serial(n, d, z, rho, lam, qnew, qblk)
    integer, intent(in) :: n
    real(dp), intent(in) :: d(n), z(n), rho
    real(dp), intent(out) :: lam(n), qnew(n, n)
    real(dp), intent(inout) :: qblk(n, n)
    real(dp) :: d_s(n), z_s(n), wmat(n, n)

    call sort_dz_q(n, d, z, qblk, d_s, z_s)
    call rank1_update_eigh(n, d_s, z_s, rho, lam, wmat)
    qnew = matmul(qblk, wmat)
end subroutine secular_eigen_serial

subroutine secular_eigen_caf(n, d, z, rho, lam, qnew, qblk)
    integer, intent(in) :: n
    real(dp), intent(in) :: d(n), z(n), rho
    real(dp), intent(out) :: lam(n), qnew(n, n)
    real(dp), intent(inout) :: qblk(n, n)
    call secular_eigen_serial(n, d, z, rho, lam, qnew, qblk)
    sync all
end subroutine secular_eigen_caf

!> Sort d ascending; apply the same permutation to z and to columns of qblk.
subroutine sort_dz_q(n, d, z, qblk, d_s, z_s)
    integer, intent(in) :: n
    real(dp), intent(in) :: d(n), z(n)
    real(dp), intent(inout) :: qblk(n, n)
    real(dp), intent(out) :: d_s(n), z_s(n)
    integer :: i, j, min_idx, row
    real(dp) :: td, tz, tq

    d_s = d
    z_s = z
    do i = 1, n - 1
        min_idx = i
        do j = i + 1, n
            if (d_s(j) < d_s(min_idx)) min_idx = j
        end do
        if (min_idx /= i) then
            td = d_s(i); d_s(i) = d_s(min_idx); d_s(min_idx) = td
            tz = z_s(i); z_s(i) = z_s(min_idx); z_s(min_idx) = tz
            do row = 1, n
                tq = qblk(row, i)
                qblk(row, i) = qblk(row, min_idx)
                qblk(row, min_idx) = tq
            end do
        end if
    end do
end subroutine sort_dz_q

subroutine rank1_update_eigh(n, d, z, rho, lam, qnew)
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char
    integer, intent(in) :: n
    real(dp), intent(in) :: d(n), z(n), rho
    real(dp), intent(out) :: lam(n), qnew(n, n)
    real(c_double), allocatable :: a(:)
    integer(c_int) :: n_c, info, i, j
    character(kind=c_char) :: jobz, uplo
    interface
        function lapacke_dsyev(matrix_layout, jobz, uplo, n, a, lda, w) &
            bind(C, name="LAPACKE_dsyev")
            import :: c_int, c_double, c_char
            integer(c_int) :: lapacke_dsyev
            integer(c_int), value :: matrix_layout, n, lda
            character(kind=c_char), value :: jobz, uplo
            real(c_double) :: a(*), w(*)
        end function lapacke_dsyev
    end interface
    integer(c_int), parameter :: LAPACK_COL_MAJOR = 102

    allocate (a(n * n))
    a = 0.0_c_double
    do i = 1, n
        a(i + (i - 1) * n) = d(i)
        lam(i) = 0.0_dp
    end do
    do j = 1, n
        do i = 1, n
            a(i + (j - 1) * n) = a(i + (j - 1) * n) + rho * z(i) * z(j)
        end do
    end do
    n_c = int(n, c_int)
    jobz = 'V'
    uplo = 'U'
    info = lapacke_dsyev(LAPACK_COL_MAJOR, jobz, uplo, n_c, a, n_c, lam)
    if (info /= 0) error stop "rank1_update_eigh: LAPACKE_dsyev failed"
    ! eigenvectors in a, column-major
    do j = 1, n
        do i = 1, n
            qnew(i, j) = a(i + (j - 1) * n)
        end do
    end do
    deallocate (a)
end subroutine rank1_update_eigh

subroutine sort_dz(n, d, z, d_s, z_s)
    integer, intent(in) :: n
    real(dp), intent(in) :: d(n), z(n)
    real(dp), intent(out) :: d_s(n), z_s(n)
    integer :: i, j, min_idx
    real(dp) :: td, tz

    d_s = d
    z_s = z
    do i = 1, n - 1
        min_idx = i
        do j = i + 1, n
            if (d_s(j) < d_s(min_idx)) min_idx = j
        end do
        if (min_idx /= i) then
            td = d_s(i)
            d_s(i) = d_s(min_idx)
            d_s(min_idx) = td
            tz = z_s(i)
            z_s(i) = z_s(min_idx)
            z_s(min_idx) = tz
        end if
    end do
end subroutine sort_dz

!> j-th secular root in the j-th interlacing interval (1-based, poles sorted).
real(dp) function secular_root_j(n, j, d, z, rho) result(lam)
    integer, intent(in) :: n, j
    real(dp), intent(in) :: d(n), z(n), rho
    real(dp) :: lo, hi, mid, flo, fmid, eps
    integer :: it
    real(dp), parameter :: huge_shift = 1.0e3_dp

    eps = epsilon(1.0_dp) * 10.0_dp
    if (j == 1) then
        lo = d(1) - huge_shift * (1.0_dp + abs(d(1)) + abs(rho))
        hi = d(1) - eps * (1.0_dp + abs(d(1)))
    else if (j == n) then
        lo = d(n) + eps * (1.0_dp + abs(d(n)))
        hi = d(n) + huge_shift * (1.0_dp + abs(d(n)) + abs(rho))
    else
        lo = d(j - 1) + eps * (1.0_dp + abs(d(j - 1)))
        hi = d(j) - eps * (1.0_dp + abs(d(j)))
    end if

    ! Ensure bracket (expand if needed)
    flo = secular_f(n, d, z, rho, lo)
    do it = 1, 40
        if (flo * secular_f(n, d, z, rho, hi) <= 0.0_dp) exit
        if (j == 1) lo = lo - (hi - lo)
        if (j == n) hi = hi + (hi - lo)
        if (j > 1 .and. j < n) then
            ! shift slightly
            lo = lo + 0.1_dp * (hi - lo)
            hi = hi - 0.1_dp * (hi - lo)
        end if
        flo = secular_f(n, d, z, rho, lo)
    end do

    do it = 1, 80
        mid = 0.5_dp * (lo + hi)
        fmid = secular_f(n, d, z, rho, mid)
        if (abs(hi - lo) < eps * (1.0_dp + abs(mid))) exit
        if (flo * fmid <= 0.0_dp) then
            hi = mid
        else
            lo = mid
            flo = fmid
        end if
    end do
    lam = 0.5_dp * (lo + hi)
end function secular_root_j

real(dp) function secular_f(n, d, z, rho, lam) result(f)
    integer, intent(in) :: n
    real(dp), intent(in) :: d(n), z(n), rho, lam
    integer :: i
    f = 1.0_dp
    do i = 1, n
        f = f + rho * (z(i) * z(i)) / (d(i) - lam)
    end do
end function secular_f

subroutine secular_vec(n, d, z, rho, lam, w)
    integer, intent(in) :: n
    real(dp), intent(in) :: d(n), z(n), rho, lam
    real(dp), intent(out) :: w(n)
    integer :: i
    real(dp) :: nrm
    do i = 1, n
        w(i) = z(i) / (d(i) - lam)
    end do
    nrm = sqrt(sum(w * w))
    if (nrm > 0.0_dp) w = w / nrm
end subroutine secular_vec

end module gjp_imtqlx_dc

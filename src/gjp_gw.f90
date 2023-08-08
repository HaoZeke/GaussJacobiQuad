module gjp_gw
use gjp_types, only: dp, gjp_sparse_matrix
use gjp_constants, only: pi
implicit none
contains

subroutine gauss_jacobi_gw(n, a, b, x, w)
    integer, intent(in) :: n
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: x(n), w(n)
    real(dp) :: zmom
    type(gjp_sparse_matrix) :: jacmat
    real(dp) :: d(n), e(n - 1), z(n, n), work(2 * n - 2)
    integer :: info, i

    interface
        subroutine DSTEQR(COMPZ, N, D, E, Z, LDZ, WORK, INFO)
            use gjp_types, only: dp
            character :: COMPZ
            integer :: N, LDZ, INFO
            real(dp) :: D(*), E(*), Z(LDZ, *), WORK(*)
        end subroutine DSTEQR
    end interface
    ! Compute the Jacobi matrix
    jacmat = jacobi_matrix(n, a, b)
    print*,"The diagonal elements"
    print*,jacmat%diagonal(1:n)
    print*,"The off-diagonal elements"
    print*,jacmat%off_diagonal(1:n - 1)
    zmom = jacobi_zeroeth_moment(a, b)

    ! Extract diagonal and off-diagonal elements
    d = jacmat%diagonal(1:n)
    e = jacmat%off_diagonal(1:n - 1)

    ! Initialize z as identity matrix
    z = 0.0_dp
    do i = 1, n
        z(i, i) = 1.0_dp
    end do

    ! Diagonalize the Jacobi matrix.
    call DSTEQR('V', n, d, e, z, n, work, info)

    if (info /= 0) then
        print*,'Error in DSTEQR:', info
        return
    end if

    ! The eigenvalues are the nodes
    x = d
    ! The weights are related to the squares of the first components of the eigenvectors
    w = z(1, :)**2 * zmom

end subroutine gauss_jacobi_gw

function jacobi_matrix(n, alpha, beta) result(jacmat)
    integer, intent(in) :: n ! Size of the matrix, number of points
    real(dp), intent(in) :: alpha, beta
    type(gjp_sparse_matrix) :: jacmat
    integer :: idx
    real(dp) :: ab, abi, a2b2

    allocate (jacmat%diagonal(n))
    allocate (jacmat%off_diagonal(n - 1))

    ab = alpha + beta
    abi = 2.0_dp + ab
    jacmat%diagonal(1) = (beta - alpha) / abi
    jacmat%off_diagonal(1) = 4.0_dp * (1.0_dp + alpha) * (1.0_dp + beta) &
                             / ((abi + 1.0_dp) * abi * abi)
    a2b2 = beta * beta - alpha * alpha
    do idx = 2, n
        abi = 2.0_dp * idx + ab
        jacmat%diagonal(idx) = a2b2 / ((abi - 2.0_dp) * abi)
        abi = abi**2
        if (idx < n) then
            jacmat%off_diagonal(idx) = 4.0_dp * idx * (idx + alpha) * (idx + beta) &
                                       * (idx + ab) / ((abi - 1.0_dp) * abi)
        end if
    end do
    jacmat%off_diagonal(1:n - 1) = sqrt(jacmat%off_diagonal(1:n - 1))
end function jacobi_matrix

function jacobi_zeroeth_moment(alpha, beta) result(zmom)
    real(dp), intent(in) :: alpha, beta
    real(dp) :: zmom
    real(dp) :: ab, abi

    ab = alpha + beta
    abi = 2.0_dp + ab

    zmom = 2.0_dp**(alpha + beta + 1.0_dp) * exp(log_gamma(alpha + 1.0_dp) &
                                                 + log_gamma(beta + 1.0_dp) - log_gamma(abi))

    if (zmom <= 0.0_dp) then
        error stop "Zeroth moment is not positive but should be"
    end if
end function jacobi_zeroeth_moment
end module gjp_gw

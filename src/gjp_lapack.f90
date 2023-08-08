module gjp_lapack
implicit none

integer, parameter :: dp = kind(0.d0)

interface
    subroutine DSTEQR(COMPZ, N, D, E, Z, LDZ, WORK, INFO)
        use gjp_types, only: dp
        character :: COMPZ
        integer :: N, LDZ, INFO
        real(dp) :: D(*), E(*), Z(LDZ, *), WORK(*)
    end subroutine DSTEQR
end interface

contains

end module gjp_lapack

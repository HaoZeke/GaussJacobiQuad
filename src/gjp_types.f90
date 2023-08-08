module gjp_types
implicit none
private
public sp, dp, hp, qp, ivector, dvector, zvector, gjp_sparse_matrix

integer, parameter :: dp = kind(0.d0), & ! double precision
                      hp = selected_real_kind(15), & ! high precision
                      qp = selected_real_kind(32), & ! quadruple precision
                      sp = kind(0.) ! single precision

type ivector ! allocatable integer vector
    integer, pointer :: vec(:) => null()
end type

type dvector ! allocatable real double precision vector
    real(dp), pointer :: vec(:) => null()
end type

type zvector ! allocatable complex double precision vector
    complex(dp), pointer :: vec(:) => null()
end type

!> @brief Sparse representation of a Jacobi matrix.
  !!
  !! @details The type jacobi_sparse_matrix represents a Jacobi matrix in a sparse format.
  !! It holds the diagonal and off-diagonal elements of the matrix, as the other elements are zeros.
  !! The sparse format is suitable for efficient storage and computation when dealing with Jacobi matrices.
type gjp_sparse_matrix
    real(dp), allocatable :: diagonal(:) !< Diagonal elements of the Jacobi matrix.
    real(dp), allocatable :: off_diagonal(:) !< Off-diagonal elements of the Jacobi matrix.
end type gjp_sparse_matrix

end module

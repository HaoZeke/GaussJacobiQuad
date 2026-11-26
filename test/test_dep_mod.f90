program test_dep_mod
  use GaussJacobiQuad, only: dp
  implicit none
  real(dp) :: x
  x = 1.0_dp
  if (x <= 0.0_dp) error stop 1
end program test_dep_mod

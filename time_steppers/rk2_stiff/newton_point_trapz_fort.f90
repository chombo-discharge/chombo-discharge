subroutine SOLVE_LU(N, J, F, INFO) bind(C, name="FORT_SOLVE_LU")
  use iso_c_binding, only : c_int, c_double
  implicit none

  integer(c_int), intent(inout)  :: N
  integer(c_int), intent(inout)  :: INFO
  real(c_double), intent(inout)  :: J(N,N)
  real(c_double), intent(inout)  :: F(N);

  integer :: i
  integer :: IPIV(N)

  ! Solving J*x = -F
  F = -1.0*F

  call dgesv(N, 1, J, N, IPIV, F, N, INFO)

end subroutine SOLVE_LU

pure subroutine newton_point_trapz(i) bind(C, name="newton_point_trapz_fort")
  use iso_c_binding, only : c_int
  implicit none
  integer(c_int), intent(inout) :: i
  i = i + 1
end subroutine newton_point_trapz

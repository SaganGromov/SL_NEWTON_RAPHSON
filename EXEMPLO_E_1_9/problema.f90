module problema
   implicit none
contains
   pure function p(x) result(val)
      real(8), intent(in) :: x
      real(8)            :: val
      val = 1.0d0
   end function p
   pure function q(x) result(val)
      real(8), intent(in) :: x
      real(8)            :: val
      val = 0
   end function q
   pure function r(x) result(val)
      real(8), intent(in) :: x
      real(8)            :: val
      val = 1
   end function r
end module problema
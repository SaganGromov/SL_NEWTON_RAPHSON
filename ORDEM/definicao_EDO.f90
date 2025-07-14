module definicao_EDO
   use problema
   implicit none
contains
   subroutine f_ode4(x,y,lam,dy)
      real(8), intent(in)  :: x,lam
      real(8), intent(in)  :: y(4)
      real(8), intent(out) :: dy(4)
      real(8) :: u,pu,w,pw
      u  = y(1);  pu = y(2)
      w  = y(3);  pw = y(4)

      dy(1) = pu / p(x)
      dy(2) = (q(x) - lam*r(x))*u
      dy(3) = pw / p(x)
      dy(4) = (q(x) - lam*r(x))*w - r(x)*u
   end subroutine f_ode4
end module definicao_EDO
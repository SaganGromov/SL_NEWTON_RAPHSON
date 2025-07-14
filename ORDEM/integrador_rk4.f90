module integrador_rk4
   use definicao_EDO
   implicit none
contains
   subroutine passo_rk4(x,y,h,lam)
      real(8), intent(inout) :: x
      real(8), intent(inout) :: y(:)
      real(8), intent(in)    :: h,lam
      integer                :: n
      real(8), allocatable   :: k1(:),k2(:),k3(:),k4(:),yt(:)

      n = size(y)
      allocate(k1(n),k2(n),k3(n),k4(n),yt(n))

      call f_ode4(x,           y,            lam, k1)
      yt = y + 0.5d0*h*k1 ; call f_ode4(x+0.5d0*h, yt, lam, k2)
      yt = y + 0.5d0*h*k2 ; call f_ode4(x+0.5d0*h, yt, lam, k3)
      yt = y +       h*k3 ; call f_ode4(x+      h, yt, lam, k4)

      y = y + h*(k1 + 2.0d0*k2 + 2.0d0*k3 + k4)/6.0d0
      x = x + h

      deallocate(k1,k2,k3,k4,yt)
   end subroutine passo_rk4
end module integrador_rk4
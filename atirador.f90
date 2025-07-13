module atirador
   use constantes
   use problema
   use integrador_rk4
   implicit none
contains
   !----------------------------------------------------------
   ! Integra de a até b e devolve:
   !   F(λ)   = γ u(b) + δ u'(b)
   !   F'(λ)  = γ w(b) + δ w'(b)
   !   nodes  = nº de zeros de u em (a,b)
   !----------------------------------------------------------
   subroutine shoot(lam, F, Fp, nodes)
      real(8), intent(in)  :: lam
      real(8), intent(out) :: F, Fp
      integer, intent(out) :: nodes
      real(8) :: y(4), x, h, yprev, v, vp, norm
      integer :: i

      ! Condições iniciais em x=a
      norm = sqrt(alpha*alpha + beta*beta)
      x    = a
      y    = 0.0d0
      y(1) = -beta  / norm
      y(2) =  alpha / norm * p(a)

      h     = (b-a)/real(passos,8)
      nodes = 0
      yprev = y(1)

      do i = 1, passos
         call passo_rk4(x,y,h,lam)
         if (y(1)*yprev < 0.0d0) nodes = nodes + 1
         yprev = y(1)
      end do

      v  = y(2)/p(b)
      vp = y(4)/p(b)

      F  = gamma*y(1) + delta*v
      Fp = gamma*y(3) + delta*vp
   end subroutine shoot
end module atirador
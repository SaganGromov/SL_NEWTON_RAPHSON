module salvador_solucoes
   use constantes
   use problema
   use integrador_rk4
   implicit none
contains


   
   !----------------------------------------------------------
   ! Calcula e salva a autofunção para o autovalor
   !----------------------------------------------------------
   subroutine salvar_autofuncao(lam)
      real(8), intent(in) :: lam
      
      real(8) :: y(4), x, h, norm_int, norm
      integer :: i, unit_file, n_pontos
      real(8), allocatable :: x_vals(:), y_vals(:)
      
      ! Número de pontos para salvar (mais denso que a integração)
      n_pontos = 1000
      allocate(x_vals(n_pontos+1), y_vals(n_pontos+1))
      
      ! Condições iniciais
      norm = sqrt(alpha*alpha + beta*beta)
      x    = a
      y    = 0.0d0
      y(1) = -beta  / norm
      y(2) =  alpha / norm * p(a)
      
      h = (b-a)/real(n_pontos,8)
      
      ! Integra e armazena valores
      x_vals(1) = x
      y_vals(1) = y(1)
      
      do i = 1, n_pontos
         call passo_rk4(x, y, h, lam)
         x_vals(i+1) = x
         y_vals(i+1) = y(1)
      end do
      
      ! Normaliza a autofunção usando integral ∫ r(x)|y|² dx
      norm_int = 0.0d0
      do i = 1, n_pontos
         norm_int = norm_int + r(0.5d0*(x_vals(i)+x_vals(i+1))) * &
                              y_vals(i)**2 * (x_vals(i+1)-x_vals(i))
      end do
      norm_int = sqrt(norm_int)
      
      if (norm_int > 0.0d0) then
         y_vals = y_vals / norm_int
      end if
      
      ! Garante que a autofunção seja positiva em x=a+h
      if (y_vals(2) < 0.0d0) then
         y_vals = -y_vals
      end if
      
      ! Salva em arquivo
      unit_file = 20
      open(unit=unit_file, file='autofuncao.dat', status='replace')
      
      write(unit_file, '(A,ES20.12)') '# Autofuncao, lambda = ', lam
      write(unit_file, '(A)') '# x            y(x)'
      write(unit_file, '(A)') '#-----------------------'
      
      do i = 1, n_pontos+1
         write(unit_file, '(2(ES20.12,2X))') x_vals(i), y_vals(i)
      end do
      
      close(unit_file)
      
      print *, 'Autofuncao salva em autofuncao.dat'
      
      deallocate(x_vals, y_vals)
   end subroutine salvar_autofuncao
   
end module salvador_solucoes

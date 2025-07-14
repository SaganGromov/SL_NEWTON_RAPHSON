module verificador_ordem
   use constantes
   use problema
   use integrador_rk4
   use atirador
   use newton_raphson
   implicit none
   
   ! For this specific problem (p=1, q=0, r=1 with y(0)=y(1)=0),
   ! the exact eigenvalues are λ_n = (nπ)² and eigenfunctions are sin(nπx)
   
contains

   !----------------------------------------------------------
   ! Computes the exact eigenfunction at point x for eigenvalue lam
   !----------------------------------------------------------
   function exact_eigenfunction(x, lam) result(y_exact)
      real(8), intent(in) :: x, lam
      real(8) :: y_exact
      real(8) :: n
      
      ! For λ = (nπ)², we have n = sqrt(λ)/π
      n = sqrt(lam) / acos(-1.0d0)
      
      ! The exact eigenfunction is sin(nπx)
      y_exact = sin(n * acos(-1.0d0) * x)/(n*acos(-1.0d0))
   end function exact_eigenfunction
   
   !----------------------------------------------------------
   ! Integrates and evaluates y at a specific point x_eval
   !----------------------------------------------------------
   recursive subroutine evaluate_at_point(lam, x_eval, n_steps, y_value)
      real(8), intent(in)  :: lam, x_eval
      integer, intent(in)  :: n_steps
      real(8), intent(out) :: y_value
      
      real(8) :: y(4), x, h, norm, y_mid
      integer :: i
      
      ! Initial conditions
      norm = sqrt(alpha*alpha + beta*beta)
      x    = a
      y    = 0.0d0
      y(1) = -beta  / norm
      y(2) =  alpha / norm * p(a)
      
      h = (b-a)/real(n_steps,8)
      
      ! Integrate up to x_eval
      do i = 1, n_steps
         if (x + h > x_eval) then
            ! Last step to reach exactly x_eval
            h = x_eval - x
            call passo_rk4(x, y, h, lam)
            exit
         end if
         call passo_rk4(x, y, h, lam)
         if (x >= x_eval) exit
      end do
      
      y_value = y(1)
      
      ! Normalize by matching at x = 0.5
      if (abs(x_eval - 0.5d0) > 1.0d-10) then
         ! Need to integrate again to get normalization
         call evaluate_at_point(lam, 0.5d0, n_steps, y_mid)
         if (abs(y_mid) > 1.0d-10) then
            y_value = y_value / abs(y_mid)
         end if
      else
         ! Already at midpoint, normalize to 1
         if (abs(y_value) > 1.0d-10) then
            y_value = y_value / abs(y_value)
         end if
      end if
   end subroutine evaluate_at_point
   
   !----------------------------------------------------------
   ! Main routine to check the order of the method
   !----------------------------------------------------------
   subroutine check_method_order()
      real(8) :: lam, lam_exact, F, Fp
      real(8) :: x_eval, y_exact, y_num, error
      integer :: nodes, i, n_h
      logical :: converged
      integer :: unit_y, unit_err
      

      integer, parameter :: n_tests = 9
      integer :: steps_array(n_tests)
      real(8) :: h_array(n_tests)
      
      ! Variables for convergence order calculation
      real(8) :: prev_error, order
      
      ! Evaluation point (choose something not at boundaries)
      x_eval = 0.2d0
      
      ! First, find the eigenvalue with very high precision
      print *, '================================================'
      print *, 'VERIFICAÇÃO DA ORDEM DO MÉTODO RK4'
      print *, '================================================'
      print *, ''
      
      ! Use the original tolerance for finding eigenvalue
      lam = lam_chute
      call refine_eigenvalue(lam, lam_exact, converged)
      
      if (.not. converged) then
         print *, 'ERRO: Não foi possível encontrar o autovalor!'
         return
      end if
      
      print *, 'Autovalor encontrado: λ =', lam_exact
      print *, 'Autovalor teórico esperado: λ₁ = π² ≈', acos(-1.0d0)**2
      print *, ''
      
      ! Exact solution at x_eval
      y_exact = exact_eigenfunction(x_eval, lam_exact)
      print *, 'Ponto de avaliação: x =', x_eval
      print *, 'Valor exato y(x) =', y_exact
      print *, ''
      
      ! Prepare step sizes (from coarse to fine)
      do i = 1, n_tests
         steps_array(i) = 10 * 2**(i-1)  ! 10, 20, 40, 80, ...
         h_array(i) = (b-a) / real(steps_array(i), 8)
      end do
      
      ! Open output files
      unit_y = 30
      unit_err = 31
      open(unit=unit_y, file='convergencia_y.dat', status='replace')
      open(unit=unit_err, file='convergencia_erro.dat', status='replace')
      
      write(unit_y, '(A)') '# Convergência de y(x) com h'
      write(unit_y, '(A,F10.6)') '# x = ', x_eval
      write(unit_y, '(A,ES20.12)') '# lambda = ', lam_exact
      write(unit_y, '(A)') '# h                    y(x)                  n_steps'
      write(unit_y, '(A)') '#----------------------------------------------------'
      
      write(unit_err, '(A)') '# Convergência do erro com h'
      write(unit_err, '(A,F10.6)') '# x = ', x_eval
      write(unit_err, '(A,ES20.12)') '# lambda = ', lam_exact
      write(unit_err, '(A,ES20.12)') '# y_exact = ', y_exact
      write(unit_err, '(A)') '# h                    |y_exact - y_num|     log2(erro_i/erro_i+1)'
      write(unit_err, '(A)') '#--------------------------------------------------------------'
      
      print *, 'Calculando para diferentes tamanhos de passo...'
      print *, ''
      print *, '    h           passos        y(x)          |erro|        ordem'
      print *, '----------------------------------------------------------------'
      
      prev_error = 0.0d0
      
      do i = 1, n_tests
         call evaluate_at_point(lam_exact, x_eval, steps_array(i), y_num)
         error = abs(y_exact - y_num)
         
         write(unit_y, '(ES20.12, 2X, ES20.12, 2X, I10)') &
               h_array(i), y_num, steps_array(i)
         
         if (i > 1 .and. prev_error > 0.0d0 .and. error > 0.0d0) then
            order = log(prev_error/error) / log(2.0d0)
            write(unit_err, '(ES20.12, 2X, ES20.12, 2X, F10.4)') &
                  h_array(i), error, order
            print '(ES12.4, I10, ES14.6, ES12.4, F10.4)', &
                  h_array(i), steps_array(i), y_num, error, order
         else
            write(unit_err, '(ES20.12, 2X, ES20.12, 2X, A)') &
                  h_array(i), error, '    ---'
            print '(ES12.4, I10, ES14.6, ES12.4, A)', &
                  h_array(i), steps_array(i), y_num, error, '      ---'
         end if
         
         prev_error = error
      end do
      
      close(unit_y)
      close(unit_err)
      
      print *, ''
      print *, 'Resultados salvos em:'
      print *, '  - convergencia_y.dat    (h vs y(x))'
      print *, '  - convergencia_erro.dat (h vs erro)'
      print *, ''
      print *, 'Para plotar no gnuplot:'
      print *, '  plot "convergencia_erro.dat" u 1:2 w lp title "Erro", \'
      print *, '       x**4 title "O(h^4)"'
      
   end subroutine check_method_order
   
end module verificador_ordem
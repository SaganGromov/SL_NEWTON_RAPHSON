module newton_raphson
   use constantes
   use atirador
   implicit none
contains
   !----------------------------------------------------------
   ! Refina um autovalor usando Newton-Raphson
   ! Entrada: lam_inicial (chute inicial)
   ! Saída:   lam_refinado, converged (se convergiu)
   !----------------------------------------------------------
   subroutine refine_eigenvalue(lam_inicial, lam_refinado, converged)
      real(8), intent(in)  :: lam_inicial
      real(8), intent(out) :: lam_refinado
      logical, intent(out) :: converged
      
      real(8) :: lam, lam_old, F, Fp
      integer :: iter, nodes, max_iter
      real(8) :: delta_lam
      
      max_iter = 100
      lam = lam_inicial
      converged = .false.
      
      do iter = 1, max_iter
         lam_old = lam
         
         ! Calcula F(λ) e F'(λ)
         call shoot(lam, F, Fp, nodes)
         
         ! Verifica se F é pequeno o suficiente
         if (abs(F) < tol_Newt) then
            converged = .true.
            lam_refinado = lam
            return
         end if
         
         ! Verifica se F' é muito pequeno (evita divisão por zero)
         if (abs(Fp) < 1.0d-14) then
            converged = .false.
            return
         end if
         
         ! Passo de Newton-Raphson
         delta_lam = -F/Fp
         lam = lam + delta_lam
         print '(A,I0,A,G0,A,G0)', 'Iteração:', iter, '  λ:', lam, '  F(λ):', F
         
         ! Verifica convergência em λ
         if (abs(delta_lam) < tol_Newt) then
            converged = .true.
            lam_refinado = lam
            return
         end if
      end do
      
      ! Se chegou aqui, não convergiu
      converged = .false.
      lam_refinado = lam
   end subroutine refine_eigenvalue
   
   !----------------------------------------------------------
   ! Procura um autovalor próximo a lam_guess
   ! Usa bissecção para encontrar mudança de sinal em F
   !----------------------------------------------------------
   subroutine find_eigenvalue_bracket(lam_guess, lam_min, lam_max, found)
      real(8), intent(in)  :: lam_guess
      real(8), intent(out) :: lam_min, lam_max
      logical, intent(out) :: found
      
      real(8) :: lam, F, Fp, F_left, F_right, delta
      integer :: nodes, i
      
      delta = 0.1d0
      found = .false.
      
      ! Procura mudança de sinal em F ao redor de lam_guess
      do i = 1, 100
         lam_min = lam_guess - i*delta
         lam_max = lam_guess + i*delta
         
         call shoot(lam_min, F_left, Fp, nodes)
         call shoot(lam_max, F_right, Fp, nodes)
         
         if (F_left * F_right < 0.0d0) then
            found = .true.
            return
         end if
         
         ! Aumenta o passo de busca
         if (i > 10) delta = 0.5d0
         if (i > 50) delta = 1.0d0
      end do
   end subroutine find_eigenvalue_bracket
   
end module newton_raphson

program sturm_liouville_solver
   use constantes
   use atirador
   use newton_raphson
   use salvador_solucoes
   implicit none
   
   real(8) :: lam, lam_refinado, F, Fp
   integer :: nodes
   logical :: converged
   
   print *, ''
   print *, 'Chute inicial: lambda =', lam_chute
   print *, ''
   print '(A,G0)', 'Tolerância para o método de Newton-Raphson definida como: ', tol_Newt
   
   ! Avalia o chute inicial
   lam = lam_chute
   call shoot(lam, F, Fp, nodes)
   
   print *, 'Avaliacao do chute inicial:'
   print '(A,G0)', '  F(lambda)  = ', F
   print '(A,G0)', '  F''(lambda) = ', Fp
   print *, ''
   
   ! Refina o autovalor usando Newton-Raphson
   print *, 'Refinando com Newton-Raphson...'
   call refine_eigenvalue(lam, lam_refinado, converged)
   
   if (converged) then
      print *, ''
      print *, 'Newton-Raphson convergiu!'
      print '(A,G0)', 'Autovalor refinado: lambda = ', lam_refinado
      
      ! Verifica o resultado
      call shoot(lam_refinado, F, Fp, nodes)
      print '(A,G0)', 'Verificacao: F(lambda) = ', F
      
      ! Salva os resultados
      print *, ''
      print *, 'Salvando resultados...'
      call salvar_autofuncao(lam_refinado)
      
      print *, ''
      print *, 'Calculo concluido com sucesso!'
   else
      print *, ''
      print *, 'ERRO: Newton-Raphson nao convergiu!'
      print *, 'Tente um chute inicial diferente.'
   end if
   
end program sturm_liouville_solver

module constantes
   implicit none
   !--------- Domínio --------------------------------------------------
   real(8), parameter :: a = 0.0d0
   real(8), parameter :: b = 1.0d0
   real(8), parameter :: lam_chute = 1.4d0  ! Chute inicial para o autovalor
   !--------- Condições de contorno  α y(a)+β y'(a)=0 , γ y(b)+δ y'(b)=0
   real(8), parameter :: alpha = 1.0d0, beta  = 0.0d0     ! y(a)=0
   real(8), parameter :: gamma = 1.0d0, delta = 0.0d0     ! y(b)=0
   !--------- Malha e controle ----------------------------------------
   integer, parameter :: passos   = 1.0d5
   real(8), parameter :: tol_Newt = 1.0d-12     
end module constantes
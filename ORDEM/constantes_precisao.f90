module constantes
   implicit none
   !--------- Domínio --------------------------------------------------
   real(8), parameter :: a = 0.0d0
   real(8), parameter :: b = 1.0d0
   real(8), parameter :: lam_chute = 1.0d1  ! Close to π² ≈ 9.8696
   !--------- Condições de contorno  α y(a)+β y'(a)=0 , γ y(b)+δ y'(b)=0
   real(8), parameter :: alpha = 1.0d0, beta  = 0.0d0     ! y(a)=0
   real(8), parameter :: gamma = 1.0d0, delta = 0.0d0     ! y(b)=0
   !--------- Malha e controle ----------------------------------------
   integer, parameter :: passos   = 1.0d6
   real(8), parameter :: tol_Newt = 1.0d-16    ! Machine precision
end module constantes

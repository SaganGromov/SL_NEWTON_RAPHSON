Para compilar e executar:
```bash
gfortran -O3 constantes_precisao.f90 problema.f90 definicao_EDO.f90 integrador_rk4.f90 atirador.f90 newton_raphson.f90 verificador_ordem.f90 main_verificacao.f90 -o verificar_ordem && ./verificar_ordem
```
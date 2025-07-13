Para compilar e executar os programas, execute o comando a seguir em suas respectivas pastas:
```bash
gfortran -O3 constantes.f90 problema.f90 definicao_EDO.f90 integrador_rk4.f90 atirador.f90 newton_raphson.f90 salvador_solucoes.f90 main.f90 -o sturm_liouville && ./sturm_liouville
```
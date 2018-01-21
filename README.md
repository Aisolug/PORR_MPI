# PORR_MPI


MPI implementacja MPICH wersja 3.0.4

Kompilacja: 

mpicc.mpich -std=c99 main.c testFun.c jacobi.c richardson.c -o program -lm

Uruchomienie:

mpiexec -np 8 ./program richardson

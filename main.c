#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include "testFun.h"
#include "jacobi.h"
#include "richardson.h"
#include "math.h"
#include <string.h>

int main(int argc, char* argv[]) {

    char *algorithm;
    algorithm = argv[1];

//zmienne potrzebne do obliczen

    funkcja pierwsza;
    jacobi jacobi;
    richardson rich;

    int n = 2024;
    double precision = 0.0000000000001;
    double precision_rich = 0.01;
    pierwsza.nr = jacobi.nr = rich.nr = 4;
    pierwsza.size = jacobi.size = rich.size = n;


/// Tworzenie funkcji testowej ///

    int *vector;

    if (pierwsza.nr == 1) {
        exponent = 3.0;
        multiplier = 3.0;
        vector = createTestF1(pierwsza);
    }

    if (pierwsza.nr == 2) {
        exponent = 5.0;
        multiplier = 5.0;
        vector = createTestF2(pierwsza);
    }

    if (pierwsza.nr == 3) {
        exponent = 4.0;
        multiplier = 4.0;
        vector = createTestF3(pierwsza);
    }

    if (pierwsza.nr == 4) {
        vector = createTestF4(pierwsza);
        exponent = 7.0;
        multiplier = 4.0;
    }
    free(pierwsza.test);

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int rank_size;
    MPI_Comm_size(MPI_COMM_WORLD, &rank_size);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    count = n/rank_size;

    double *abs = (double *) calloc(count, sizeof(double));
    double *prev = (double *) calloc(jacobi.size, sizeof(double));

    if (argc < 2) {
        if (rank == 0) {
            printf("\n No algorithm name!\n");
            printf("\n for Jacobi write: jacobi\n");
            printf(" for Richardson write: richardson\n\n");
        }
    } else if (strcmp("jacobi", algorithm) == 0) {

        double *jacobian;
        jacobi.fun = (double *) calloc(jacobi.size, sizeof(double));

        int msec_jacobi;
        clock_t start_jacobi, diff_jacobi;
        if (rank == 0) start_jacobi = clock();

        for (loop = 1; 2 < 3; loop++) {

            int sum = 0;
            int out = 0;

            jacobian = resultJacobi(jacobi, vector, jacobian);

            if (loop > 1) {
                int i;
                int iter = 0;
                for (i = rank * count; iter < count && out == 0; i++, iter++) {
                    abs[iter] = fabs(jacobian[i] - prev[i]);
                    if (abs[iter] > precision) out = 1;
                }
                MPI_Allreduce(&out, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                if (sum == 0) {
                    free(jacobi.fun);
                    break;
                }
                    //printf("Wykonane iteracje dla Jacobiego: %d, rank %d\n", loop, rank);
            }

            int p;
            for (p = 0; p < jacobi.size; p++) {
                prev[p] = jacobian[p];
            }
        }
        if (rank == 0) {
            diff_jacobi = clock() - start_jacobi;
            msec_jacobi = diff_jacobi * 1000 / CLOCKS_PER_SEC;
            printf("\nCzas wykonania dla Jacobiego: %d\n", msec_jacobi % 1000);
        }
    }
    else if (strcmp("richardson", algorithm) == 0) {

        double* richardson;
        rich.fun = (double *) calloc(rich.size, sizeof(double));

        int msec_richardson;
        clock_t start_richardson, diff_richardson;
	    if (rank == 0) start_richardson = clock();

        for (loop = 1; 2 < 3; loop++){

            int sum = 0;
            int out = 0;

	        richardson = resultRichardson(rich, vector, richardson);

            if (loop > 1) {
                int i;
                int iter = 0;
                for (i = rank * count; iter < count && out == 0; i++, iter++) {
                    abs[iter] = fabs(richardson[i] - prev[i]);
                    if (abs[iter] > precision_rich) out = 1;
                }
                MPI_Allreduce(&out, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
                if (sum == 0) {
                    free(rich.fun);
                    break;
                }
            }

            int p;
            for (p = 0; p < rich.size; p++) {
                prev[p] = richardson[p];
            }

        }
        if (rank == 0) {
            diff_richardson = clock() - start_richardson;
            msec_richardson = diff_richardson * 1000 / CLOCKS_PER_SEC;
            printf("\nCzas wykonania dla Richardsona: %d\n", msec_richardson % 1000);
        }
    }
    else{
        if (rank == 0) {
            printf("\n Wrong algorithm name!\n");
            printf("\n for Jacobi write: jacobi\n");
            printf(" for Richardson write: richardson\n\n");
        }
    }

    free(abs);
    free(prev);

    MPI_Finalize();
}

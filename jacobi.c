#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "jacobi.h"
#include <mpi.h>

double* resultJacobi(jacobi jacobi, int* test, double* jacobian){

    double* in = (double *) calloc(count, sizeof(double));
    jacobi.input = (double *) calloc(jacobi.size, sizeof(double));

	if (loop == 1){

		int z;
        int iter = 0;
		for (z = rank*count  ; iter < count; z++,iter++){
			in[iter] = 2.0;
		}
    }

	else{
        int z;
        int iter = 0;
        for (z = rank*count  ; iter < count; z++,iter++){
            in[iter] = jacobian[z];
        }
	}

    MPI_Gather(in, count, MPI_DOUBLE, jacobi.input, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(jacobi.input, jacobi.size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int i;
    int iter1=0;
    jacobi.sum = (double *) calloc(jacobi.size, sizeof(double));
    double* sumP = (double *) calloc(count, sizeof(double));

    for (i = rank*count; iter1 < count; i++, iter1++){
        int k;
		for (k = 0; k < jacobi.size; k ++){
            if(i != k){
                sumP[iter1] -= (double)(test[i*jacobi.size+k])*jacobi.input[k];
            }
        }
    }

    MPI_Gather(sumP, count, MPI_DOUBLE, jacobi.sum, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(jacobi.sum, jacobi.size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* funP = (double *) calloc(count, sizeof(double));

	int n;
    int iter = 0;
    for (n = rank*count  ; iter < count; n++,iter++){
        funP[iter] = pow((jacobi.sum[n] - (n+1)*multiplier)/(test[n+n*jacobi.size]*jacobi.input[n]),1.0/exponent);
    }

    MPI_Gather(funP, count, MPI_DOUBLE, jacobi.fun, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(jacobi.fun, jacobi.size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(in);
    free(funP);
    free(sumP);
    free(jacobi.input);
    free(jacobi.sum);

	return jacobi.fun;
}



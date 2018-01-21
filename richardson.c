#include <stdio.h>
#include <stdlib.h>
#include "richardson.h"
#include <math.h>

double* resultRichardson(richardson rich, int* test, double* richardson){

	double gamma = 0.0000003;
    double* in = (double *) calloc(count, sizeof(double));
	rich.input = (double *) calloc(rich.size, sizeof(double));
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
            in[iter] = richardson[z];
        }
	}

    MPI_Gather(in, count, MPI_DOUBLE, rich.input, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(rich.input, rich.size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    rich.sum = (double *) calloc(rich.size, sizeof(double));
    double* sumP = (double *) calloc(count, sizeof(double));

    int i;
    int iter1 = 0;
    for (i = rank; iter1 < count; i++, iter1++){
	int k;
		for (k = 0; k < rich.size; k ++){
            if(i != k){
                sumP[iter1] += (double)(test[i*rich.size+k])*rich.input[k];
            }
        }
    }

	MPI_Gather(sumP, count, MPI_DOUBLE, rich.sum, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(rich.sum, rich.size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double* funP = (double *) calloc(count, sizeof(double));
	int n;
	int iter = 0;

    for (n = rank*count; iter < count; n++,iter++){
        funP[iter] = rich.sum[n] + (double)test[n+n*rich.size]*pow(rich.input[n],exponent) + (n+1)*multiplier;
        funP[iter] = rich.input[n] - gamma*funP[iter];
    }
	MPI_Gather(funP, count, MPI_DOUBLE, rich.fun, count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(rich.fun, rich.size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(sumP);
	free(funP);
	free(rich.input);	
	free(rich.sum);
	return rich.fun;
}


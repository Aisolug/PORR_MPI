#ifndef JACOBI_H
#define JACOBI_H
#include <mpi.h>


typedef struct tag_jacobi{
	
	int nr, size;
	double *input; 
	double *sum;
	double *fun;

}jacobi;

double* resultJacobi(jacobi jacobi, int* test, double* jacobian);


int count, loop, rank;
double exponent, multiplier;
#endif

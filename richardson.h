#ifndef RICHARDSON_H
#define RICHARDSON_H
#include <mpi.h>

typedef struct tag_richardson{
	
	int nr, size;
	double *input; 
	double *sum;
	double *fun;

}richardson;

double* resultRichardson(richardson rich, int* test, double* richardson);

int count, loop, rank;
double exponent, multiplier;
#endif

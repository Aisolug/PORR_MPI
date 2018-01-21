#include <stdio.h>
#include <stdlib.h>
#include "testFun.h"

int x;

int* complementMatrix(funkcja fun){
    	
	int i;    
        for (i = 1; i < fun.size; i ++){
		int k;            
		for (k = 0; k < fun.size; k ++){
	        	if(k < i){
	                    fun.test[i*fun.size+k] = fun.test[k*fun.size + i];    
                	}
            	}
        }
	return fun.test;
}

int* derivativeMatrix(funkcja fun, int exponent){
                
		
	int i;
        for (i = 0; i < fun.size; i++){
            fun.test[i*(1+fun.size)] = fun.test[i*(1+fun.size)]*exponent;
        }

	return fun.test;
}
	
int* createTestF1(funkcja fun){
	
	int exponent = 4;
	fun.test = (int *) calloc(fun.size*fun.size, sizeof(int));
	int coefficiency[5] = {7,20,5,-9,-10};
	x = 3;
	int j;

	for (j = x-1; j< fun.size; j++){
		int i;
            	for (i = 1; i < x; i++){
                	fun.test[(j-i)*fun.size + j-i] += coefficiency[x-1-i];
           	}
            	fun.test[j+j*fun.size] += j*coefficiency[2];
            	
		for (int i = 0; i < x-1; i++){
                	fun.test[(j-2)*fun.size +j-i] = (j+1)*coefficiency[x*2-2-i];
 	        }  
        }
	fun.test[0] += 5; //zapewniona dominacja
	fun.test = complementMatrix(fun);
	fun.test = derivativeMatrix(fun, exponent);
	return fun.test;
}

int* createTestF2(funkcja fun){
	
	int exponent = 6;
	fun.test = (int *) calloc(fun.size*fun.size, sizeof(int));
	int coefficiency[5] = {4,5,2,-6,-7};
	x = 3;
	int j;

	for (j = x-1; j< fun.size; j++){
		int i;
            	for (i = 1; i < x; i++){
                	fun.test[(j-i)*fun.size + j-i] += coefficiency[x-1-i];
           	}
            	fun.test[j+j*fun.size] += (j+1)*coefficiency[2];
            	
		for (int i = 0; i < x-1; i++){
                	fun.test[(j-2)*fun.size +j-i] = j*coefficiency[x*2-2-i];
 	        }  
        }

	fun.test = complementMatrix(fun);
	fun.test = derivativeMatrix(fun, exponent);
	return fun.test;
}
     
int* createTestF3(funkcja fun){
	
	int exponent = 5;
	fun.test = (int *) calloc(fun.size*fun.size, sizeof(int));
	int coefficiency[7] = {5,15,8,3,-5,-6,-7};
	x = 4;
	int j;

	for (j = x-1; j< fun.size; j++){
		int i;
            	for (i = 1; i < x; i++){
                	fun.test[(j-i)*fun.size + j-i] += coefficiency[x-1-i];
           	}
            	fun.test[j+j*fun.size] += (j+1)*coefficiency[3];
            	
		for (int i = 0; i < x-1; i++){
                	fun.test[(j-3)*fun.size +j-i] = j*coefficiency[x*2-2-i];
 	        }  
        }

	fun.test = complementMatrix(fun);
	fun.test = derivativeMatrix(fun, exponent);
	return fun.test;
}


int* createTestF4(funkcja fun){
	
	int exponent = 8;
	fun.test = (int *) calloc(fun.size*fun.size, sizeof(int));
	int coefficiency[7] = {20,25,30,5,-15,-16,-17};
	x = 4;
	int j;

	for (j = x-1; j< fun.size; j++){
		int i;
            	for (i = 1; i < x; i++){
                	fun.test[(j-i)*fun.size + j-i] += coefficiency[x-1-i];
           	}
            	fun.test[j+j*fun.size] += (j+1)*coefficiency[3];
            	
		for (int i = 0; i < x-1; i++){
                	fun.test[(j-3)*fun.size +j-i] = j*coefficiency[x*2-2-i];
 	        }  
        }

	fun.test = complementMatrix(fun);
	fun.test = derivativeMatrix(fun, exponent);
	return fun.test;
}


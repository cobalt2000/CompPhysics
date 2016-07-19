#include <stdlib.h>
#include <stdio.h>

#include "eigenVec.h"


/**************************************************
This is a small program to test the library 
eigenVec.  It will create a small matrix, and
pass it to the library to solve the eigen problem.
**************************************************/
int main(int argc, char** argv){

int i,j, k;  // Standard indexs
int info;    // This is a flag used to tell if the function worked correctly.

int n = 4;
double *A; // The matrix to the solved.
double *eVectors; // Pointer to the solution. 
double *eValues;  // Pointer to the eigen values.


/* allocate memory for the matrix and the eigenvectors */
A = (double*)malloc(sizeof(double) * n * n);
eVectors = (double*)malloc(sizeof(double) * n * n);
eValues = (double*)malloc(sizeof(double) * n);


// Fill in matrix

for( i = 0; i < n; i++){
	for(j = i; j < n; j++){
		A[i + j * n] = i + j + 1;
		A[j + i * n] = i + j + 1;
	}
}
printVector(A, n * n);



// Call solver
 info = find_eigen_vectors( A, n, eValues, eVectors );
if(info != 0){
printf("There was an error computing the eigenvalues: %d\n", info);
}

// Print solution
printf("The eigen values are:\n");
printVector(eValues, n);


printf("The eigen vectors are the columns of:\n");
printVector(eVectors, n * n);



// free memory

free(A);
free(eVectors);
free(eValues);


return 0;

}


#ifndef EIGENVEC
#define EIGENVEC

int find_eigen_vectors(double **A, int n, double *eVal, double **eVec );
/*
        int find_eigen_vectors(double **A, int n, double *eVal, double ** eVec)
This funtion finds the eigen values and eigen vectors of the symmetric real matrix A.
Funtion will return 0 if no error.

Arguments:
        double **A      input           Symmetric real matrix
        int n           input           size of A
        double *eVal    output          eigen values returned in a array of size n
        double **eVec   output          egien vectors returned in the columns of the n x n array
*/


void printArray(double **A, int n, int m); // Prints the n x m array

void printVector(double *v, int n);	// Prints the vector of size n

/* Some complilers require a profile for fortran functions.
   This may be required on some computers.
*/
/*
extern "C" void dsyev( char* jobz, 
		   char* uplo, 
		   int* n, 
		   double* a, 
		   int* lda, 
		   double* w, 
		   double* work, 
		   int* lwork, 
		   int* info );
*/
#endif




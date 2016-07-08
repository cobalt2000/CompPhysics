#include <stdlib.h>
#include <stdio.h>
#include "eigenVec.h"


//#define PRINT_ALL


int find_eigen_vectors(double **A, int n, double *eVal, double **eVec ){
/* 
	int find_eigen_vectors()
This funtion finds the eigen values and eigen vectors of the symmetric real matrix A.
Funtion will return 0 if no error.

Arguments:
	double **A	input		Symmetric real matrix
	int n		input		size of A
	double *eVal	output		eigen values returned in a array of size n
	double **eVec	output		egien vectors returned in the columns of the n x n array
*/

	int i, j;
	double *U; // This is a temp matrix to pass data to LAPACK.
	int lwork = -1; 	// the size of the workspace needed by LAPACK.	
	double wkopt;
	double* work;	//pointer to workspace needed by LAPACK
	int info;

	U = (double*)malloc(sizeof(double)*n*n);
	for(i = 0; i < n; i++){
		U[i + n*i ] = A[i][i];
		for( j = (i+1); j < n; j++){
			U[i + j*n ] = A[i][j];
			U[i*n + j ] = 0.0;
		}
	}

#ifdef PRINT_ALL
	printf("The matrix passed to LAPACK is:\n");
	printVector(U, n*n);

#endif

	/*find the optimum size of workspace */
	dsyev_( "V", "U", &n, U, &n, eVal, &wkopt, &lwork, &info );
	
	lwork = (int)wkopt;
	work = (double*)malloc(sizeof(double)*lwork);

	/* Solve eigen problem */
	dsyev_( "V", "U", &n, U, &n, eVal, work, &lwork, &info );
/*
From http://www.math.utah.edu/software/lapack/lapack-d/dsyev.html
 SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK,
                        INFO )

          CHARACTER     JOBZ, UPLO

          INTEGER       INFO, LDA, LWORK, N

          DOUBLE        PRECISION A( LDA, * ), W( * ), WORK( * )

 PURPOSE
      DSYEV computes all eigenvalues and, optionally, eigenvectors
      of a real symmetric matrix A.

 ARGUMENTS
      JOBZ    (input) CHARACTER*1
              = 'N':  Compute eigenvalues only;
              = 'V':  Compute eigenvalues and eigenvectors.

      UPLO    (input) CHARACTER*1
              = 'U':  Upper triangle of A is stored;
              = 'L':  Lower triangle of A is stored.

      N       (input) INTEGER
              The order of the matrix A.  N >= 0.

      A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
              On entry, the symmetric matrix A.  If UPLO = 'U',
              the leading N-by-N upper triangular part of A con-
              tains the upper triangular part of the matrix A.  If
              UPLO = 'L', the leading N-by-N lower triangular part
              of A contains the lower triangular part of the
              matrix A.  On exit, if JOBZ = 'V', then if INFO = 0,
              A contains the orthonormal eigenvectors of the
              matrix A.  If JOBZ = 'N', then on exit the lower
              triangle (if UPLO='L') or the upper triangle (if
              UPLO='U') of A, including the diagonal, is des-
              troyed.

      LDA     (input) INTEGER
              The leading dimension of the array A.  LDA >=
              max(1,N).

      W       (output) DOUBLE PRECISION array, dimension (N)
              If INFO = 0, the eigenvalues in ascending order.

      WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)

              On exit, if INFO = 0, WORK(1) returns the optimal
              LWORK.

      LWORK   (input) INTEGER
              The length of the array WORK.  LWORK >= max(1,3*N-
              1).  For optimal efficiency, LWORK >= (NB+2)*N,
              where NB is the blocksize for DSYTRD returned by
              ILAENV.

      INFO    (output) INTEGER
              = 0:  successful exit
              < 0:  if INFO = -i, the i-th argument had an illegal
              value
              > 0:  if INFO = i, the algorithm failed to converge;
              i off-diagonal elements of an intermediate tridiago-
              nal form did not converge to zero.
*/


/* Now the problem should be solved.  The data needs to be placed in the outgoing arrays */

for( i = 0; i < n; i++){
	for(j = 0; j < n; j++){
		eVec[i][j] = U[i + j*n ];
	}
}



free(U);
free(work);

return info;

}

////////////////////////////////////////////////////////////

void printArray(double **A, int n, int m){
        int i, j;
        for( i = 0; i < n; i++){
                for(j = 0; j < m; j++){
                        printf("%1.3lg ", A[i][j]);
                }
                printf("\n");
        }

        return;
}

/////////////////////////////////////////////////////////////
void printVector(double *v, int n){
	int i;
	for(i = 0; i < n; i++)
		printf(" %2.5lg\n", v[i]);
	return;
}


#include <iostream>
#include <stdio.h>
#include <mpi.h>		/* MPI header file                         */

#include "sprng_cpp.h"
#include "eigenVec.h"


#define SEED 985456376


int main(int argc, char *argv[])
{
  int streamnum, nstreams;
  Sprng *stream;
  double rn;
  int i, myid, nprocs, j, k;
  int n = 3;
  int gtype;  /*---    */
  
  double *A, *eVec, *eVal;
  
  /*************************** MPI calls ***********************************/

  MPI_Init(&argc, &argv);	/* Initialize MPI                          */
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);	/* find process id                 */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs); /* find number of processes      */

  /************************** Initialization *******************************/

  streamnum = myid;	
  nstreams = nprocs;		/* one stream per processor                */
/*--- node 0 is reading in a generator type */
  if(myid == 0)
  {
      if(argc < 2){
          gtype = 1;
      }else{
          gtype = atoi(argv[1]);
    }
      
}
  MPI_Bcast(&gtype,1,MPI_INT,0,MPI_COMM_WORLD );
    
  stream = SelectType(gtype);
  stream->init_sprng(streamnum,nstreams,SEED,SPRNG_DEFAULT);	/* initialize stream */
  
  A = (double*)malloc(sizeof(double) * n * n);
  eVec = (double*)malloc(sizeof(double) * n * n);
  eVal = (double*)malloc(sizeof(double) * n);
  
  //printf("\n\nProcess %d, print information about stream:\n", myid);
  //stream->print_sprng();

  /*********************** print random numbers ****************************/

//  for (i=0;i<3;i++)
//  {
//    rn = stream->sprng();		/* generate double precision random number */
//    printf("Process %d, random number %d: %.14f\n", myid, i+1, rn);
//  }

  /*************************** free memory *********************************/

  /******************* Fill Random Matrix *******************************/
  for( i = 0; i < n; i++){
        A[i + i*n] = stream->sprng();
        for(j = 0; j < i; j++){
            A[i + j*n] = stream->sprng();
            A[i*n + j] = A[i + j*n];
        }
      
  }
  
  find_eigen_vectors(A, n, eVal, eVec );
  
  
  stream->free_sprng();		/* free memory used to store stream state  */

  free(A);
  free(eVec);
  free(eVal);
  
  
  MPI_Finalize();		/* Terminate MPI                           */

  return 0;
}

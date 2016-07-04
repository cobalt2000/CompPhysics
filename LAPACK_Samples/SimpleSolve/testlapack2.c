#include <stdio.h>
#include <stdlib.h>

//#include <atlas_enum.h>
//#include "clapack.h"

int main()
{
  int *ipiv;
  int i, j;
  int info;
int n = 3;
int columnsOfx = 1;

double *m;
double *x;

m = (double*)malloc(sizeof(double)*3*3);
x = (double*)malloc(sizeof(double)*3);
ipiv = (int*)malloc(sizeof(int)*3);


// Store by col, not by row.  That is, this matrix is labled as:
//   [ 0 3 6 ] 
//   [ 1 4 7 ]
//   [ 2 5 8 ]
//
// 


m[0] = 3.0;
m[1] = 1.0;
m[2] = 2.0;
m[3] = 1.0;
m[4] = 5.0;
m[5] = 6.0;
m[6] = 3.0;
m[7] = 9.0;
m[8] = 5.0;

x[0] = -1.0;
x[1] = 3.0;
x[2] = -3.0;


  for (i=0; i<3; ++i) {
    for (j=0; j<3; ++j)  printf("%5.1f", m[i+j*3]);
    putchar('\n');
  }
  dgesv_(&n, &columnsOfx, m, &n, ipiv, x,&n , &info );
  if (info != 0) fprintf(stderr, "failure with error %d\n", info);

  for (i=0; i<3; ++i) printf("%5.1f %3d\n", x[i], ipiv[i]);


free(m);
free(x);
free(ipiv);

  return 0;
}

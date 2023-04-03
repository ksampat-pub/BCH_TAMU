/*******************************************************
This routine solves a system of simultaneous equations of any order using the direct methdo called Gaussian inversion.
Note that the power representations are kept intact even when solving linear equations unlike in the case of lu_decomp.c where pow to poly conversion takes place.  
********************************************************/
#include<stdio.h>
#include"global.h"
extern void gau_lin_solve(int MM[][T],int degr,int *X,int *B)
{
  int TT;
  int i,j,k,sum;
  void rearrange(int a[][T],int b[T], int uu); 
  TT = degr;
  /* COnvert to lower triangular form  */
  for( i=((int)T-TT);i<T;i++) {
    if( MM[i][i] == -1)
      rearrange(MM,B,i);
    for( j=i+1;j<T;j++) {
      if( MM[j][i]!= -1) {
	for( k=i+1;k<T;k++) {
	  MM[j][k] = gf_minus(gf_mult(gf_mult(MM[j][k],inv[ MM[j][i] ]),MM[i][i]),MM[i][k]);
	}
	B[j] = gf_minus(gf_mult(gf_mult( B[j],inv[ MM[j][i] ]),MM[i][i]),B[i]);
      }
    }
  }
  /* Solve */
  for(i=((int)T-1);i>=((int)T-TT);i--) {
    sum = -1; 
    for(k=((int)T-1);k>i;k--) 
      sum=gf_add(sum,gf_mult(X[k],MM[i][k]));
    X[i] =gf_mult( gf_minus( B[i],sum),inv[ MM[i][i] ] );
  }
}
void rearrange( int MM[][T], int *B , int col)
{
  int i,temp,j;
  for(i=col+1;i<T;i++) {
    if( MM[i][col] != -1) {
      for(j=col;j<T;j++) {
	temp = MM[col][j];
	MM[col][j] = MM[i][j];
	MM[i][j] = temp;
      }
      temp = B[col];
      B[col] = B[i];
      B[i] = temp;
      break;
    }
  }
}


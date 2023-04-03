/* Routine to calculate determinant of square matrix of order T  */

#include<stdio.h>
#include"global.h"
/* int *parent ;*/            /*parent is global  */
int parent[T];
extern void determ( int MM[][T],int TT)
{
  extern int det_sum;
  int det ( int , int bs[][T], int,int);
  if( (TT == 1) && ( MM[T-1][T-1] == -1)) {
      det_sum = N+1;            /* An illegal value for det_sum */
      return;
  }
  if( TT == 1)
    return;
  det ( 0 ,MM , TT,1);
}

int det (int curr_row,int  MM[][T],int  TT, int first_time)
{
  extern int det_sum;
  static int ii,temp_sum;
  int i,k,jj,kk,a,b,c,d,flag=0;
  i=0;
  if ( first_time ) {
    det_sum = temp_sum = 0;
    ii = 1;
    for(i=0;i<T;i++)
      *(parent+i) = 0;
  }

  /*penultimate and ultimate row operation */
  if ( (T-curr_row) == 2) {
    for(k=T-TT;k<T;k++) {
      if ( parent[k] != 1 ) {
	if( flag == 0) {
	  flag = 1;
	  a = MM[TT-2][k];
	  c = MM[TT-1][k];
        }
	else {
	  flag = 0;
	  b = MM[TT-2][k];
	  d = MM[TT-1][k];
	}  
      }
    }
    if(TT==T && first_time ) det_sum=gf_minus( gf_mult(a,d),gf_mult(b,c));
    return gf_minus( gf_mult(a,d),gf_mult(b,c));
  }

  /* main recursive function that calculates the determinant */
  for(jj=T-TT;jj<T;jj++)    {
    if ( jj == T-TT) kk = 0;
    if ( parent[jj] != 1) {
      ii = pow( (-1),kk++);
      parent[jj] = 1;
      curr_row = curr_row+1;
      if(ii==-1)
	temp_sum = gf_minus(temp_sum , gf_mult(det(curr_row,MM,TT,0),MM[curr_row-1][jj]));
      if(ii==1)
	temp_sum = gf_add(temp_sum , gf_mult(det(curr_row,MM,TT,0),MM[curr_row-1][jj])); 
      parent[jj] = 0;
      curr_row = curr_row - 1 ;
      if(curr_row == 0 ) {
	det_sum = gf_add(det_sum , temp_sum);
	temp_sum = -1;  /**/
      }
    }
  }
  return temp_sum ;
}







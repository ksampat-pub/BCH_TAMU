#include<stdio.h>
#include"global.h"
int gf_add(int a, int b)
{
  if(a==-1) return b;
  if(b==-1) return a;
  if(a==0) a=N;
  if(b==0) b=N;
  return Add_tab[a][b];
}

int gf_mult ( int a, int b)
{
  if( (a == -1 )||( b== -1 )) return -1;
  else return (a+b)%N;    /* mult implies addition of powers  */
}

int gf_pow( int a, int b)
{
  if (a==-1) return -1;
/*   if (b==-1) return 1; */
  if(b==0) return 0;
  return ((a*b) % N) ; /* returns the power representation of the element */ 
}

int gf_minus(int a, int b)
{
  int i;
  if( (a== -1) && (b== -1) ) 
    return -1;
  if( b== -1 ) 
    return a;
  if( a==b ) 
    return -1;
  if( b == 0 ) 
    b =N;
  for(i=1;i<= N;i++)
    if(Add_tab[b][i] == a ) {
      if(i==N) 
	i=0;
      return i;
    }
}
/* Multiplies the message polynomial with the generator polynomial*/
void gf_poly_mult(int *code, int dc, int *m, int dm)
  {
 int temp[N][N];
 int i,j;
 for(i=0;i<dm;i++) {
   for(j=0;j<dc;j++) {
     if(j<i) 
       temp[i][j] = -1;
     if( (i+j) < N ){
       if(j<N-K+1)
	 temp[i][i+j] = gf_mult(m[i],gen[j]);
       if( j>=N-K+1 )
	 temp[i][i+j] = -1;
     }
   }
 }
 for(j=0;j<dc;j++) { 
   code[j] = -1;
   for(i=0;i<dm;i++) {
     code[j]= gf_add(code[j],temp[i][j]);
   }
 }
}

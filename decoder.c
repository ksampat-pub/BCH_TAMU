#include"global.h"
#include<stdio.h>
int det_sum;
int decode ( int ch2[N] , int ch1[N] )
{
  /*prog to decode a BCH  code over GF(q^S) */
  int k,deg,X[T],E[T],cnt=0, Lambda[2][2*T+1], L, Tx[2][2*T+1], Delta,FLAG1,FLAG2;
  int i,j,sndr[2*T+1],sum,S1[T],S2[T],B[T][T],A[T][T],LL[T];
  int gau_lin_solve(int MM[][T],int orde,int *a,int *b);
  void determ( int PP[][T],int ord);
  det_sum=-1;
 
  /*r(x) is of length N in prog we go from lsb to msb*/
 
  
   i=0; 
   for(j=0;j<N;j++)       { 
     if(ch1[j]!=ch2[j]) 
       i++;                  /* Simplifies debugging by counting errors */ 
   } 
   if( i==0)       {   
     return 0;   
   } 
   sndr[0] = -1;
   for(j=1;j<=2*T;j++)       {
     sum=-1;
     for(i=0;i<N;i++)
       sum = gf_add( sum,gf_mult( ch1[i],gf_pow(element[j+R],i) ) );/*modd*/
     sndr[j] = sum;
     if( sndr[j]==-1)
       cnt++;
     /*calculates the syndromes*/
   }
   if(cnt==2*T) {
     /* check if all syndromes are zero*/
     /*This means that the original CW has been transformed */
     /* into another codeword. Thf. undetectable error */
     return -1;
   }
   if ( PGZ == 0) {
     /* Berlekamp Massey Decoder        */
     
     k = 0;
     L = 0;
     for ( i=0;i<=2*T;i++) {
       Tx[0][i] = -1;
       Tx[1][i] = -1;
       Lambda[0][i] = -1;
       Lambda[1][i] = -1;
     }
     Lambda[0][0] = 0;
     Tx[0][1] = 0;
     FLAG1 = FLAG2 = 0;
     while ( k < 2*T ) {
       sum = -1;
       k++;
       for ( i = 1; i<=L ; i++) 
	 sum = gf_add(sum, gf_mult(Lambda[0][i],sndr[k-i]));
       
       Delta = gf_minus( sndr[k], sum);
       if ( Delta == -1 )
	 goto Eight;
       if ( FLAG1&&FLAG2 || !(FLAG1&&FLAG2) ) {
	 for ( i = 0;i <=k;i++)
	   Lambda[1][i] = gf_minus( Lambda[0][i] , gf_mult(Delta,Tx[0][i]));
	 FLAG1=0;
       }
       if ( 2*L >= k)
	 goto Eight;
       FLAG1 = 1 ;
       L = k - L;
       for ( i = 0;i<k;i++) {
	 Tx[0][i+1] =gf_mult(Lambda[0][i], (int)(N-Delta)%N );
       }
       
     Eight:	  FLAG2 = 1 ;
     for ( i = k+1;i>=1;i--) {
       Tx[1][i] = Tx[0][i-1];
       if ( Delta != -1)
	 Lambda[0][i] = Lambda[1][i]; 
       Tx[1][0] = -1;
     }
     
     if ( !FLAG1 || Delta == -1 ) 
       for ( i=0; i<=k; i++) 
	 Tx[0][i] = Tx[1][i];
     if (  Delta != -1 )
       Lambda[0][0] = Lambda[1][0]; 
     
     }
     
     k=0;
     for(i=1;i<=N;i++) {
       sum = 0;
       for(j=1;j<= T;j++) {
	 sum = gf_add( sum, gf_mult(Lambda[0][j],gf_pow(element[i],j)));
       }
       if(sum==-1) {
	 X[k] = element[i];  /* Inverse of error locator*/
	 /* Solution to eqn. 9-20 Wicker  */
	 k++;
       }	
     }
     /* End of Berlekamp Massey decoder */
   }
   
   
   else {
     /* The Peterson Gorenstein Zierler decoder  */ 
     /* Determinant et al.   */
     for(i=1;i<=T;i++) {
       for(j=0;j<T;j++) 
	 A[i-1][j] = sndr[i+j];
       S1[i-1] =  gf_minus(-1, sndr[i+T] );
     }
     deg = T;
     while( (det_sum==-1) && (deg>0)) { 
       determ(A,deg); 
       if( (deg == 1) && (det_sum == N+1) )
	 deg = -2;
       deg-- ;
     }
     deg++;
     if(deg < 0)
       return -1;
     gau_lin_solve( A, deg,LL,S1);
     k=0;
     for(i=1;i<=N;i++) {
       sum = 0;
       for(j=1;j<=deg;j++) {
	 sum = gf_add( sum, gf_mult(LL[T-j],gf_pow(element[i],j)));
       }
       if(sum==-1) {
	 X[k] = element[i];  /* Inverse of error locator*/
	 /* Solution to eqn. 9-20 Wicker  */
	 k++;
       }	
     }  
   }
   
   if(k==0) 
     return -1;
   
   for(i=0;i<k;i++)
     for(j=0;j<k;j++)
       if(i!=j)
	 if(X[i]==X[j]) 
	   return -1;
   /* if roots not found or repeated roots found, then uncorrectable error */   
   
   /*Now to calculate the error locations  */
   /* First we form the B matrix ( B of eq 9-26) */
   
   for(i=1;i<=k;i++) {
     for(j=0;j<k;j++) {
       B[i-1+T-k][j+T-k] = gf_pow((N-X[j])%N,i+R);
     }
     S2[i-1+T-k] = sndr[i];   /* T-deg term to make data cons. w/ det o/p */
   }
   gau_lin_solve(B,k,E,S2);
   
   for(i=0;i<k;i++) 
     ch1[locator[X[i]]] = gf_minus(ch1[locator[X[i]]],E[T-k+i]);
   
   /* Final Check  */
   for(j=0;j<N;j++)        {
     if(ch1[j]!=ch2[j])
       {
	 return -1;}
   }
   return 1;
   
}      /*  <------*/










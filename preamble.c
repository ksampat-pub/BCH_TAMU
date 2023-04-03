/*********************************************************************
This file generates the length of the codeword, the number of info bits required from the info given in the defines. It also generates the addition table, the generator polynomial given the minimal polunomial. The global variables are defined in the file global.h 
***************************************************************/
#include<stdio.h>
#include<math.h>
#include"global.h"
 void preamble ( void )
{

  int i,j,k,m,cnt,temp[N+1],select,primitive=0,*monic;
  int PRIM_POLY[]={1,0,0,0,0,1,1};
  /*Store in the form A0 + A1X + A2X^2 + A3X^3 .... */ 

  element = (int *)malloc( (N+1)*sizeof(int));  /* remember to free  */ 
  locator = (int *)malloc( N*sizeof(int));  /* remember to free  */
  gen = (int *)malloc( (N-K+1)*sizeof(int));/* rem to free       */
  monic = (int *)malloc( (2*S*T)*sizeof(int)); /* remember to free  */
 
  element[0]= -1;
  for(i=1;i<N;i++)    
    element[i]  = i;
  element[N] = 0;
  /* Elements represented as powers of alpha, the Nth root of unity */
  /* Notation: alpha^-1=0, alpha^0=alpha^N=1                        */

  if ( (M==1) && (S==1) ) {
    for(i=1;i<Q;i++) {
      for(j=1;j<Q;j++) {  
	if( (int)pow(i,j)%(int)Q == 1) {
	  if(j==(Q-1)) { 
	    primitive = i;
	    break;
	  }
	  else break;
	}
      }
    }
    
    
    /* code below arranges elements in order of power of the primitive */
    /* element which is more amenable for standardisation              */
    for(i=1;i<Q;i++) 
      elem[i][0] = (int)pow(primitive,i)%(int)Q;
    elem[0][0] = 0;
    
    /* Following segment generates the addition table  */
    for (i=1;i<=N;i++) {
      for(j=1;j<=N;j++) {
	temp[0] = ( elem[i][0]+elem[j][0] )%Q;
	for(m=0;m<=N;m++) 
	  if(temp[0] == elem[m][0]) 
	    select = m;
	if(select == 0) select--;
	if(select == N) select = 0;
	Add_tab[i][j] = select;
      }
    }
    
  }
  else {
    /* Elements are stored as powers of the primitive element */
    /* When S not equal to one we let beta be the primitive element */
    /*   that satisfies the primitive polynomial, order of beta is taken to*/
    /*be Q-1 but we store 0 at elem[Q-1],(in the hope for consistency*/
    /* Thf "element" stores the power representation of each element*/
    /* We use "elem" to store the poly rep. of each element        */
    /* This is done to facilitate generation of addition table     */
    
    for (i=1;i<=N;i++) {
      for(j=0;j<DEG;j++) {
	if ( i<DEG) {
	  elem[0][j] = 0;
	  elem[i][j] = 0;
	  elem[i][i] = 1;
	}
	if ( i == DEG )  { 
	  elem[i][j]=((int)P-PRIM_POLY[j])%(int)P;
	}
        if ( i > DEG ) {
	  if( j >=1) 
	    elem[i][j] =((elem[i-1][DEG-1]*(P-PRIM_POLY[j]) )%P+ elem[i-1][j-1] )%P;
	  if ( j==0)
	    elem[i][j] = (elem[i-1][DEG-1]*(P-PRIM_POLY[0]))%P;
	}
      }
    }
    
    /* Following segment generates the addition table  */
    for (i=1;i<=N;i++) {
      for(j=1;j<=N;j++) {
	for(k=0;k<DEG;k++)
	  temp[k] = ( elem[i][k]+elem[j][k] )%P;
	for(m=0;m<=N;m++) {
	  for(k=0;k<DEG;k++) {
	    if(temp[k] == elem[m][k]) 
	      if( k== (DEG - 1) ) 
		select = m;
	      else continue;
	    else break;
	  }
	}
	if(select == 0) select--;
	if(select == N) select = 0;
	Add_tab[i][j] = select;
      }
    }
  }
  /*Because of the simple representation the locator is simpler  */
  for(i=0;i<N;i++) 
    locator[i] = (N-i)%N;
  
  /* Forming the generator after forming cyclotomic cosets  */
  
  monic[0] = 1+R;
  j=1;
  k=1;
  while ( ((int)(element[1+R]*pow(Q,k))%(int)N) != element[1+R]) {
    for( m=0; m<j; m++ ) {
      if( monic[m]!=monic[j] ) {
	if ( m == j-1 ) {
	  monic[j] = (int)(element[1+R]*pow(Q,k))%(int)N;
	  k++;
	}
      }
      else
	break;
    }
    j++;
  }
  for( i=2;i<=2*T;i++) { 
    k = 1;
    monic[j] = i+R;
    for( m=0; m<j; m++ ) {
      if( monic[m]!=monic[j] ) {
	if ( m == j-1 ) {
	  monic[j] = i+R;
	}
      }
      else
	break;
    }
    if ( m!= j ) 
      continue;
    j++;
    while ( ((int)(element[i+R]*pow(Q,k))%(int)N) != element[i+R]) {
      for( m=0; m<j; m++ ) {
	if( monic[m]!=monic[j] ) {
	  if ( m == j-1 ) {
	    monic[j] = (int)(element[i+R]*pow(Q,k))%(int)N;
	    k++;
	  }
	}
        else
	  break;
      }
      j++;
    }
  }
  K = N-j;  


  /* generating the generator! (polynomial) */ 
  for ( cnt=0; cnt <= N-K; cnt++ )  
    gen[cnt] = -1; 
  gen[0] = 0;                    /* stands for unity  */
  for ( cnt=1; cnt <= N-K; cnt++ ) { 
    for(i=N-K; i>0 ; i--)  
      gen[i] = gf_add( gen[i-1], gf_mult( gf_minus(-1,monic[cnt-1]) , gen[i] )); 
    gen[0] =                     gf_mult( gf_minus(-1,monic[cnt-1]) , gen[0] ); 
    
  } 
  
  /* gen[0]  =  27;
  gen[1]  =  45;
  gen[2]  =   0;
  gen[3]  =  -1;
  gen[4]  =  -1;
  gen[5]  =  27;
  gen[6]  =  18;
  gen[7]  =  36;
  gen[8]  =   0; */
  
  /* To form the inversion matrix over GF(Q)  */
  
  /*   for(i=1;i<Q;i++)  */
  /*     for(j=1;j<Q;j++)  */
  /*     if( (i*j)%Q == 1 ) */
  /*       inv[i] = j; */
  for(i=1;i<N;i++)
    inv[i] = N-i;
  inv[0] = 0; 
  free(monic);
}  /* Preamble ends here.  */


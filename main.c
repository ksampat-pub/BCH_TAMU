/********************************************************************* 
SEE GLOBAL.H BEFORE MODIFYING. 
Main routine to evaluate the performance of BCH codes over GF(q^s) 
where q,s are any integers with certain restrictions in AWGN channel
************************************************************************/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<sys/types.h>
#include<unistd.h>
#include"global.h"

main(int argc, char *argv[])
{
  extern int *element,*gen;
  int *message,i,j,k,code[N],final[N],cnt_error=0,index,nearest[DEG],l,m,n;
  int low_lim, up_lim,recovered[DEG];
  float mult,no_of_elems,co_ord[2],new_co_ord[2],u,v,b,y1,y2,min_dist,dist;
  unsigned long cnt,decimal_sum;
  double prob_error, SNR=0.0, Eb,Es,err_mag ;
  FILE *fp;

  /* Constellation */
  float cordx[]= {0,1,0.6235,-0.2225,-0.9010,-0.9010,-0.2225, 0.6235} ;
  float cordy[]= {0,0,0.7818, 0.9749, 0.4339,-0.4339,-0.9749,-0.7818} ;
  /* The above two lines will be required to be modified when  parameters*/
  /* in the global.h file change  */

  char *prog = argv[0];            /* name of the program for errors  */
  message = (int *)malloc( K*sizeof(int));  /* rem to free this  */
  preamble();
  Es = 7.0/8.0;                      /*depends on constellation  */
  Eb = Es*N*log10(2)/(K*log10(Q));               
 /*  srand48(getpid());  */           /* Commented for Testing */
  /* Error check */
  if ( K > N ) {
    fprintf ( stderr, "%s: K is greater than N \n", prog);
    exit(2);
  }


  mult = 10;
  no_of_elems = N+1;
  while( no_of_elems/(float)10 > 1) {
    no_of_elems =no_of_elems/10;
    mult = mult*10;
  }
  while(SNR<14) {
    cnt = 0;
    err_mag = sqrt ( Eb/( 2*pow(10, (SNR/10) ) ) );
    while ( cnt_error < 20 ) {
      i=0;
      while(i<K) {
	while( (index =(int)mult*drand48()) > N  ); 
	if( index%(int)((pow(Q,S)-1)/(Q-1)) == 0) {                    /* code is over GF(q) */
	  message[i] = element[index];
	  i++;
	}
      }
      
      gf_poly_mult(code,N,message,K);
    
      /* do mapping from codeword to element here */
       if(DIM>1) {
	   for(j=0;j<N;j++)       {
	     low_lim = 0; 
	     up_lim = (DEG/DIM);
	     while( up_lim <= DEG ) {
	       decimal_sum = 0;
	       for(i=low_lim;i<up_lim;i++)       {
		 k = i-low_lim+1;
		 decimal_sum = decimal_sum + elem[code[j]+1][i]*pow(P,(DEG/DIM)-k);
	       }
	       co_ord[0]= cordx [ decimal_sum ]; /*x coord*/
	       co_ord[1]= cordy [ decimal_sum ]; /*y co ord*/
	       u = drand48();
	       v = drand48();
	       b=sqrt(-2*log(u));
	       y1=err_mag*b*cos(2*3.14159*v); /* indep. Gaussian ns.*/
 	       y2=err_mag*b*sin(2*3.14159*v); /* y1 and y2 ~N(0,s^2)*/
	       /* add noise to mapped element here */
	       new_co_ord[0]=co_ord[0]+y1;
	       new_co_ord[1]=co_ord[1]+y2;
	       min_dist=10*(err_mag+1);
	       /* recover received  word */	
	       for( l=0; l<pow( P,(DEG/DIM ) ) ;l++ )        {
		 dist=  sqrt(pow((new_co_ord[0]-cordx[l]),2)
			     +pow((new_co_ord[1]-cordy[l]),2));
		 if(dist<min_dist)       {
		   min_dist = dist;
		   nearest[0]=l;
		 }
	       }
	        for( l=0;l<(DEG/(float)DIM-1) ;l++ )        {
		  if(( nearest[0] - pow(P, DEG/(float)DIM -l - 1)) > 0)
		    nearest[0]=recovered[up_lim-l-1] = nearest[0]%(int)pow(P, DEG/(float)DIM -l - 1);
		  else   
		    recovered[up_lim-l-1] = 0;
		}
		recovered[(int)(up_lim-(DEG/(float)DIM))] = nearest[0];
		low_lim+=(DEG/DIM)  ;
		up_lim+=(DEG/DIM);
	     }
	   
	     for(m=0;m<=N;m++) {
	       for(n=0;n<DEG;n++) {
		 if( elem[m][n]== recovered[n]) 
		   if(n==(DEG-1))
		     final[j]=m-1;
		   else;
		 else
		   break;
	       }
	     }
	   }
       }
       
       else {
	 for(j=0;j<N;j++)       {
	   k=0;
	   if( code[j] == -1) {
	     co_ord[0]= cordx [ 0  ]; /*x coord*/
	     co_ord[1]= cordy [ 0  ]; /*y co ord*/
	   }
	   else if( code[j] == 0 ) {
	     co_ord[0]= cordx [ Q-1 ];
	     co_ord[1]= cordy [ Q-1 ];
	   }
	   else {
	     co_ord[0]= cordx [  code[j]/(int)((pow(Q,S)-1)/(Q-1)) ];  /* modd area*/
	     co_ord[1]= cordy [ code[j]/(int)((pow(Q,S)-1)/(Q-1))  ];
	   }		      
	   u = drand48();
	   v = drand48();
	   b = sqrt(-2*log(u));
	   y1 = err_mag*b*cos(2*3.14159*v); /* indep. Gaussian ns.*/
	   y2 = err_mag*b*sin(2*3.14159*v);/* y1 and y2 ~N(0,s^2)*/
	   /* add noise to mapped element here */
	   new_co_ord[0]=co_ord[0]+y1;
	   new_co_ord[1]=co_ord[1]+y2;
	   min_dist=10*(err_mag+1);
	   /* recover received  word */	
	   for(l=0;l<Q ;l++)       {                /* modd area*/
	     dist=  sqrt(pow((new_co_ord[0]-cordx[l]),2)
			 +pow((new_co_ord[1]-cordy[l]),2));
	     if(dist<min_dist)       {
	       min_dist = dist;
	       nearest[k]=l;               /*modd*/
	     }
	   }
	   if( nearest[k] == 0) 
	     final[j] = -1;      
	   else if( nearest[k] == Q-1 )
	     final[j] = 0;
	   /*modd-modd means modified */
	   else 
	     final[j]= nearest[k]*(pow(Q,S)-1)/(Q-1) ;
	 }
       }
       /* demodulate */    if( cnt_error==0 && cnt > 10000000) goto niklo;
       /*       final[0]=-1   ;   For testing purpose 
       final[1]= -1  ;
       final[2]=  5 ;
       final[3]= 0  ;
       final[4]= 2  ;
       final[5]= -1  ;
       final[6]= 2  ;       */
       n=decode(code,final);
       if(n<0) 
	 cnt_error++;
       cnt++;
    }
    prob_error=(double)cnt_error/(float)cnt;
    fp = fopen("BCHQ8S2T15R2.dat","a+");
    /* Open appropriate file */
    if ( SNR == 0 )
      fprintf( fp, "%c (%d,%d,%d) (N,K,T) code \n",'%',N,K,T);
    fprintf(fp,  "%11.9f  %11.9f \n",SNR,prob_error);
    fclose (fp);
    cnt_error=0;
    SNR = SNR + 1;
  }
niklo:  free(message);
  free(gen);
  free(element);
  free(locator);

}

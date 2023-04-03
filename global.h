/******************************************************************
Changes required to specify different code are:
Name        File
P           global.h
M           global.h
Q           global.h
S           global.h
T           global.h
N           global.h
K           global.h
DEG         global.h
DIM         global.h
PGZ         global.h
PRIM_POLY   preamble.c 
cordx[]     main.c
cordy[]     main.c 
Es          main.c         Energy per symbol depending on co-ordinates x&y
Finally change the output of the Makefile.
********************************************************************/

#include<stdio.h>
#define P 2
#define M 3  
#define Q 8    /* Q = P^M*/
#define S 2    /* N = Q^S -1 for a type of BCH code */
#define T 15   /* Errors that can be corrected  */ 
#define N 63   /* Q^S -1  */
/* K is as unknown before the formation of gen poly */
#define DIM 1
#define DEG 6  /* degree of PRIM_POLY= M*S   */
#define R 2    /* To generate non primitive BCH codes */
#define PGZ 0  /* PGZ = 0 chooses the Berlekamp Massey algo. PGZ = 1 chooses the PGZ algo*/

/*Store in the form A0 + A1X + A2X^2 + A3X^3 ... in preamble.c */
int *gen,Add_tab[N+1][N+1],*element,elem[N+1][DEG],*locator,inv[N+1],K;
int gf_mult(int a, int b), gf_add(int a, int b), gf_pow(int a, int b); 
int gf_minus(int a, int b); 
void gf_poly_mult( int *w, int z,int *x,int y),preamble(void); 

/*User Defined*/
int decode( int *A,int *B); 


 
















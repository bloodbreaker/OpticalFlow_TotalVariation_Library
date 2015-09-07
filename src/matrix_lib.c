/*****************************************************************************/
/*                                                                           */
/*     Copyright 08/2006 by Dr. Andres Bruhn, Prof. Joachim Weickert,        */
/*                       and Oliver Demetz                                   */
/*     Faculty of Mathematics and Computer Science, Saarland University,     */
/*                           Saarbruecken, Germany.                          */
/*                                                                           */
/* This version has been assembled for MVTec in accordance with the          */
/* conditions specified in the contract. No further guarantees are given     */
/* beyond the ones in the contract.                                          */
/*                                                                           */
/*                                                                           */
/*****************************************************************************/


#ifndef MATRIX_LIB_INCLUDED
#define MATRIX_LIB_INCLUDED

#include <math.h>
#include "matrix_analysis_lib.c"

/* ------------------------------------------------------------------------- */

void set_matrix_2d

(
                    /*********************************************************/
    float**  A,     /* in+out : matrix (changed)                             */
    int      nx,    /* in     : matrix size in x direction                   */
    int      ny,    /* in     : matrix size in y direction                   */
    int      bx,    /* in     : boundary in x direction                      */
    int      by,    /* in     : boundary in y direction                      */
    float    a      /* in     : new value                                    */
                    /*********************************************************/
)
/* set all entries of a 2-D array A to a given value a */
{
                    /*********************************************************/
    int  i,j;       /* loop variables                                        */
                    /*********************************************************/
    

    for (i=bx; i<nx+bx; i++)
	for (j=by; j<ny+by; j++)
	{
	    A[i][j] = a;
	}
    
    return;
} 


/* ------------------------------------------------------------------------- */
void star_mul_matrix_2d

(
                    /*********************************************************/
    float**  A,     /* in	  : matrix                              */
    float**  B,     /* in	  : matrix */
    float**	 C,		/* out	  : output matrix*/
    int      nx,    /* in     : matrix size in x direction                   */
    int      ny,    /* in     : matrix size in y direction                   */
    int      bx,    /* in     : boundary in x direction                      */
    int      by    /* in     : boundary in y direction                      */
                    /*********************************************************/
)
/* set all entries of a 2-D array A to a given value a */
{
                    /*********************************************************/
    int  i,j;       /* loop variables                                        */
                    /*********************************************************/


    for (i=bx; i<nx+bx; i++)
	for (j=by; j<ny+by; j++)
	{
	    C[i][j] = A[i][j]*B[i][j];
	}

    return;
}
/* ------------------------------------------------------------------------- */

void set_matrix_twin_2d

(
                    /*********************************************************/
    float**  A1,    /* in+out : matrix (changed)                             */
    float**  A2,    /* in+out : matrix (changed)                             */
    int      nx,    /* in     : matrix size in x direction                   */
    int      ny,    /* in     : matrix size in y direction                   */
    int      bx,    /* in     : boundary in x direction                      */
    int      by,    /* in     : boundary in y direction                      */
    float    a      /* in     : new value                                    */
                    /*********************************************************/
)
/* set all entries of two 2-D array A1 and A2 to a given value a */
{
                    /*********************************************************/
    int  i,j;       /* loop variables                                        */
                    /*********************************************************/

    
    for (i=bx; i<nx+bx; i++)
	for (j=by; j<ny+by; j++)
	{
	    A1[i][j] = a;
	    A2[i][j] = a;
	}
    
    return;
} 


/* ------------------------------------------------------------------------- */

void sub_matrix_2d

(
                    /*********************************************************/
    float**  A1,    /* in     : matrix A_1                                   */
    float**  A2,    /* in     : matrix A_2                                   */
    float**  A,     /* out    : matrix A = A_1 - A_2                         */
    int      nx,    /* in     : matrix size in x direction                   */
    int      ny,    /* in     : matrix size in y direction                   */
    int      bx,    /* in     : boundary in x direction                      */
    int      by     /* in     : boundary in y direction                      */
                    /*********************************************************/
 )
/* substract the 2-D array A_2 elementwise from the 2-D array A_1 */
/* the solution is then written in the 2-D array A                */
{
                    /*********************************************************/
    int  i,j;       /* loop variables                                        */
                    /*********************************************************/


    for (i=bx; i<nx+bx; i++)
	for (j=by; j<ny+by; j++)	  
	{
	    A[i][j] = A1[i][j]-A2[i][j];
	}

    return;
} 

/* ------------------------------------------------------------------------- */

void sub_matrix_twin_2d

(
                    /*********************************************************/
    float**  A1,    /* in     : matrix A_1                                   */
    float**  B1,    /* in     : matrix B_1                                   */
    float**  A2,    /* in     : matrix A_2                                   */
    float**  B2,    /* in     : matrix B_2                                   */
    float**  A,     /* out    : matrix A = A_1 - A_2                         */
    float**  B,     /* out    : matrix B = B_1 - B_2                         */
    int      nx,    /* in     : matrix size in x direction                   */
    int      ny,    /* in     : matrix size in y direction                   */
    int      bx,    /* in     : boundary in x direction                      */
    int      by     /* in     : boundary in y direction                      */
                    /*********************************************************/
)

/* substract the 2-D array A_2 elementwise from the 2-D array A_1 */
/* the solution is then written in the 2-D array A                */
/* substract the 2-D array B_2 elementwise from the 2-D array B_1 */
/* the solution is then written in the 2-D array B                */

{
  
                    /*********************************************************/
    int  i,j;       /* loop variables                                        */
                    /*********************************************************/


    for (i=bx; i<nx+bx; i++)
	for (j=by; j<ny+by; j++)	  
	{
	    A[i][j] = A1[i][j]-A2[i][j];
	    B[i][j] = B1[i][j]-B2[i][j];
	}

    return;
}


/* ------------------------------------------------------------------------- */

void add_matrix_2d

(
                    /*********************************************************/
    float**  A1,    /* in     : matrix A_1                                   */
    float**  A2,    /* in     : matrix A_2                                   */
    float**  A,     /* out    : matrix A = A_1 + A_2                         */
    int      nx,    /* in     : matrix size in x direction                   */
    int      ny,    /* in     : matrix size in y direction                   */
    int      bx,    /* in     : boundary in x direction                      */
    int      by     /* in     : boundary in y direction                      */
                    /*********************************************************/
)

/* add the 2-D array A_1 elementwise to the 2-D array A2 */
/* the solution is then written in the 2-D array A       */

{
                    /*********************************************************/
    int  i,j;       /* loop variables                                        */
                    /*********************************************************/


    for (i=bx; i<nx+bx; i++)
	for (j=by; j<ny+by; j++)	  
	{
	    A[i][j] = A1[i][j]+A2[i][j];
	}
    
    return;
} 


/* ------------------------------------------------------------------------- */

void add_matrix_twin_2d

(
                    /*********************************************************/
    float**  A1,    /* in     : matrix A_1                                   */
    float**  B1,    /* in     : matrix B_1                                   */
    float**  A2,    /* in     : matrix A_2                                   */
    float**  B2,    /* in     : matrix B_2                                   */
    float**  A,     /* out    : matrix A = A_1 + A_2                         */
    float**  B,     /* out    : matrix B = B_1 + B_2                         */
    int      nx,    /* in     : matrix size in x direction                   */
    int      ny,    /* in     : matrix size in y direction                   */
    int      bx,    /* in     : boundary in x direction                      */
    int      by     /* in     : boundary in y direction                      */
                    /*********************************************************/
)

/* add the 2-D array A_1 elementwise to the 2-D array A2 */
/* the solution is then written in the 2-D array A       */
/* add the 2-D array B_1 elementwise to the 2-D array B2 */
/* the solution is then written in the 2-D array B       */

{
                    /*********************************************************/
    int  i,j;       /* loop variables                                        */
                    /*********************************************************/

    
    for (i=bx; i<nx+bx; i++)
	for (j=by; j<ny+by; j++)	  
	{
	    A[i][j] = A1[i][j]+A2[i][j];
	    B[i][j] = B1[i][j]+B2[i][j];
	}

    return;
} 

/* ------------------------------------------------------------------------- */

void swap_matrix_2d   

(
                     /********************************************************/
    float**  A, /* in     : matrix                                      */
    float**  B, /* out    : matrix                                      */
    int      nx,     /* in     : matrix size in x direction                  */
    int      ny,     /* in     : matrix size in y direction                  */
    int      bx,     /* in     : boundary in x direction                     */
    int      by      /* in     : boundary in y direction                     */
                     /********************************************************/
)

/* copy the matrix A_orig elementwise into the matrix A_copy */

{
                    /*********************************************************/
    int  i,j;       /* loop variables                                        */
                    /*********************************************************/
	float swap;


    for (i=bx; i<nx+bx; i++)
	for (j=by; j<ny+by; j++)	  
	{
	    swap = A[i][j];
		A[i][j] = B[i][j];
		B[i][j] = swap;
	}
    
    return;
} 

 

/* ------------------------------------------------------------------------- */


void copy_matrix_2d   

(
                     /********************************************************/
    float**  A_orig, /* in     : matrix                                      */
    float**  A_copy, /* out    : matrix                                      */
    int      nx,     /* in     : matrix size in x direction                  */
    int      ny,     /* in     : matrix size in y direction                  */
    int      bx,     /* in     : boundary in x direction                     */
    int      by      /* in     : boundary in y direction                     */
                     /********************************************************/
)

/* copy the matrix A_orig elementwise into the matrix A_copy */

{
                    /*********************************************************/
    int  i,j;       /* loop variables                                        */
                    /*********************************************************/


    for (i=bx; i<nx+bx; i++)
	for (j=by; j<ny+by; j++)	  
	{
	    A_copy[i][j] = A_orig[i][j];
	}
    
    return;
} 

 

/* ------------------------------------------------------------------------- */
void copy_matrix_twin_2d   

(
                      /*******************************************************/
    float**  A_orig,  /* in     : matrix                                     */
    float**  B_orig,  /* in     : matrix                                     */
    float**  A_copy,  /* out    : matrix                                     */
    float**  B_copy,  /* out    : matrix                                     */
    int      nx,      /* in     : matrix size in x direction                 */
    int      ny,      /* in     : matrix size in y direction                 */
    int      bx,      /* in     : boundary in x direction                    */
    int      by       /* in     : boundary in y direction                    */
                      /*******************************************************/
)

/* copy the matrix A_orig elementwise into the matrix A_copy */
/* copy the matrix B_orig elementwise into the matrix B_copy */

{
                    /*********************************************************/
    int  i,j;       /* loop variables                                        */
                    /*********************************************************/

    
    for (i=bx; i<nx+bx; i++)
	for (j=by; j<ny+by; j++)	  
	{
	    A_copy[i][j] = A_orig[i][j];
	    B_copy[i][j] = B_orig[i][j];
	}
    
    return;
} 


/* ------------------------------------------------------------------------- */

void copy_matrix_var_2d   

(
                     /********************************************************/
    float**  A_orig, /* in     : matrix                                      */
    float**  A_copy, /* out    : matrix                                      */
    int      nx,     /* in     : matrix size in x direction                  */
    int      ny,     /* in     : matrix size in y direction                  */
    int      bx1,    /* in     : boundary in x direction of A_orig           */
    int      by1,    /* in     : boundary in y direction of A_orig           */
    int      bx2,    /* in     : boundary in x direction of A_copy           */
    int      by2     /* in     : boundary in y direction of A_copy           */
                     /********************************************************/
)

/* copy the 2-D array A_orig with boundaries of size bx1 and by1 into the */
/* 2-D array A_copy with boundaries of size bx2 and by2                   */

{
                    /*********************************************************/
    int  i,j;       /* loop variables                                        */
                    /*********************************************************/


    for (i=bx1; i<nx+bx1; i++)
	for (j=by1; j<ny+by1; j++)	  
	{
	    A_copy[i-bx1+bx2][j-by1+by2] = A_orig[i][j];
	}

    return;
} 

/* ------------------------------------------------------------------------- */

void fill_matrix_pattern_2d   

(
                    /*********************************************************/
    float**  A,     /* in     : matrix                                       */
    int      nx,    /* in     : matrix size in x direction                   */
    int      ny,    /* in     : matrix size in y direction                   */
    int      bx,    /* in     : boundary in x direction                      */
    int      by     /* in     : boundary in y direction                      */
                    /*********************************************************/
)

/* fill the 2-D array A with a linear pattern */

{
 
                    /*********************************************************/
    int  i,j;       /* loop variables                                        */
    int  s,k;       /* parameters                                            */
                    /*********************************************************/
  

    /* define pattern parameters */
    k=6;    
    s=255/(k-1);
      
    /* fill array */
    for (i=bx; i<nx+bx; i++)
	for (j=by; j<ny+by; j++)	  
	{
	    A[i][j] = (i-bx)/(nx-1.0)*127+(j-by)/(ny-1.0)*127;
	}
    
    return;
} 



/* ------------------------------------------------------------------------- */

int check_matrix_2d

(
                    /*********************************************************/
    float**  A,     /* in+out : matrix                                       */
    int      nx,    /* in     : matrix size in x direction                   */
    int      ny,    /* in     : matrix size in y direction                   */
    int      bx,    /* in     : boundary in x direction                      */
    int      by     /* in     : boundary in y direction                      */
                    /*********************************************************/
)

/* check 2-D array A elementwise for nan or inf values */
{

                    /*********************************************************/
    int  i,j;       /* loop variables                                        */
                    /*********************************************************/


    for (i=bx; i<nx+bx; i++)
	for (j=by; j<ny+by; j++)
	{
	    if(isfinite(A[i][j]) == 0) return -1;
	}

    return 0;
} 

/* ------------------------------------------------------------------------- */

int scale_matrix_2d
(
		/*********************************************************/
		float**  A1,     /* in     : matrix                      */
		float**  A2,     /*    out : matrix                      */
		int      nx,     /* in     : matrix size in x direction  */
		int      ny,     /* in     : matrix size in y direction  */
		int      bx,     /* in     : boundary in x direction     */
		int      by,     /* in     : boundary in y direction     */
		float    s       /* in     : factor                      */
		/*********************************************************/
)

		/* check 2-D array A elementwise for nan or inf values */
{
int  i,j;

	for (i=bx; i<nx+bx; i++)
		for (j=by; j<ny+by; j++)
	{
		A2[i][j] = s*A1[i][j];
	}

	return 0;
} 

/* ------------------------------------------------------------------------- */

int add_scalar_matrix_2d
(
		/*********************************************************/
		float**  A1,     /* in     : matrix                      */
		float**  A2,     /*    out : matrix                      */
		int      nx,     /* in     : matrix size in x direction  */
		int      ny,     /* in     : matrix size in y direction  */
		int      bx,     /* in     : boundary in x direction     */
		int      by,     /* in     : boundary in y direction     */
		float    s       /* in     : number                      */
		/*********************************************************/
)
{
	int  i,j;

	for (i=bx; i<nx+bx; i++)
		for (j=by; j<ny+by; j++)
	{
		A2[i][j] = s+A1[i][j];
	}

	return 0;
} 

/* ------------------------------------------------------------------------- */

void abs_matrix_2d
		(
						/*********************************************************/
		float**  A,     /* in     : matrix                                       */
		float**  B,     /* out    : matrix                                       */
		int      nx,    /* in     : matrix size in x direction                   */
		int      ny,    /* in     : matrix size in y direction                   */
		int      bx,    /* in     : boundary in x direction                      */
		int      by     /* in     : boundary in y direction                      */
						/*********************************************************/
		)
{
	int i,j;
	for (i=bx; i<nx+bx; i++)
		for (j=by; j<ny+by; j++)
		{
			B[i][j] = fabsf(A[i][j]);
		}
}			

/* ------------------------------------------------------------------------- */

void rescale_affine_matrix_2d
(
						/*********************************************************/
		float**  A,     /* in     : matrix                                       */
		float**  B,     /* out    : matrix                                       */
		int      nx,    /* in     : matrix size in x direction                   */
		int      ny,    /* in     : matrix size in y direction                   */
		int      bx,    /* in     : boundary in x direction                      */
		int      by     /* in     : boundary in y direction                      */
						/*********************************************************/
)
   // rescales A to the range [0,255]
{
	int i,j;
	float minval,maxval,factor;
	minval=min_val(A,nx,ny,bx,by);
	maxval=max_val(A,nx,ny,bx,by);
	factor = 255.0f / (maxval-minval) ;
	for (i=bx; i<nx+bx; i++)
		for (j=by; j<ny+by; j++)
		{
			B[i][j] = (A[i][j]-minval)*factor;
		}
}

#endif

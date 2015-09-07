/*****************************************************************************/
/*                                                                           */
/*       Copyright 03/2012 by Dr. Andres Bruhn, Oliver Demetz and Yan Zhang  */
/*     Faculty of Mathematics and Computer Science, Saarland University,     */
/*                           Saarbruecken, Germany.                          */
/*                                                                           */
/*****************************************************************************/



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "malloc_lib.c"
#include "bounds_lib.c"
#include "funct_lib.c"
#include "matrix_lib.c"
#include "mg_trans_lib.c"
#include "diffusivity_lib.c"
#include "horn_schunck_warp_c.h"

/* ------------------------------------------------------------------------- */

void gradint_central (
    float **u,
    float **ux,
    float **uy,
    int nx,
    int ny,
    int bx,
    int by,
    float hx,
    float hy )
{
  int i,j;
  float hx_1 = 1.0/(2.0*hx);
  float hy_1 = 1.0/(2.0*hy);
  for (i = bx; i < bx+nx; i++)
    for (j = by; j < by+ny; j++){
      ux[i][j] = ( u[i+1][j] - u[i-1][j])*hx_1;
      uy[i][j] = ( u[i][j+1] - u[i][j-1])*hy_1;
      
    }
      
}



void AllocateMem2D(
                   float ***img,      // input image
                   int nx,             // x-size
                   int ny )            // y-size

{
  int i;
  
  *img = (float **) malloc(nx * sizeof(float *));
  
  for (i=0; i<nx; i++)
    (*img)[i] = (float *) malloc(ny * sizeof(float));
  
}


void DisallocateMem2D(
                      float **img,         // input memory block
                      int nx,               // x-size
                      int ny )              // y-size

{
  int i;
  
  for (i=0; i<nx; i++)
    free(img[i]);
  
  free(img);
}



/* ------------------------------------------------------------------------- */

void backward_registration
(
                           /**************************************************/
    float **f1,            /* in  : 1st image                                */
    float **f2,            /* in  : 2nd image                                */
    float **f2_bw,         /* out : 2nd image (motion compensated)           */
    float **u,             /* in     : x-component of displacement field     */
    float **v,             /* in     : y-component of displacement field     */
    int   nx,              /* in     : size in x-direction                   */
    int   ny,              /* in     : size in y-direction                   */
    int   bx,              /* in     : boundary size in x-direction          */
    int   by,              /* in     : boundary size in y-direction          */
    float hx,              /* in     : grid spacing in x-direction           */
    float hy               /* in     : grid spacing in y-direction           */
                           /**************************************************/
    )

/* creates warped version of image f2 by means of bilinear interpolation */

{
                           /**************************************************/
int   i,j;                 /* loop variables                                 */
int   ii,jj;               /* pixel coordinates                              */
float ii_fp,jj_fp;         /* subpixel coordinates                           */
float delta_i,delta_j;     /* subpixel displacement                          */
float hx_1,hy_1;           /* time saver                                     */
                           /**************************************************/

/* compute time savers */
hx_1=1.0/hx;
hy_1=1.0/hy;

/* set boundaries zero */ 
set_bounds_2d(f2,nx,ny,bx,by,(float)0.0);

 for (i=bx; i<nx+bx; i++)
     for (j=by; j<ny+by; j++)
     {		
	 /* compute subpixel location */
	 ii_fp=i+(u[i][j]*hx_1);
	 jj_fp=j+(v[i][j]*hy_1);
	               
	 /* if the required image information is out of bounds */
	 if ((ii_fp<bx)||(jj_fp<by)||(ii_fp>(nx+bx-1))||(jj_fp>(ny+by-1)))
	   {
	       /* assume zero flow, i.e. set warped 2nd image to 1st image */
	       f2_bw[i][j]=f1[i][j];	       	       	       
	      
	   }
	 /* if required image information is available */
	 else
	   {
	       /* compute integer index of upper left pixel */
	       ii=(int)floor(ii_fp);
	       jj=(int)floor(jj_fp);
	     
	       /* compute subpixel displacement */
	       delta_i = ii_fp-(float)ii;
	       delta_j = jj_fp-(float)jj;
	     
	       /* perform bilinear interpolation */
	       f2_bw[i][j]   = (1.0-delta_i)*(1.0-delta_j) * f2[ii  ][jj  ]
		                 +      delta_i *(1.0-delta_j) * f2[ii+1][jj  ]
		                 + (1.0-delta_i)*     delta_j  * f2[ii  ][jj+1]
		                 +      delta_i *     delta_j  * f2[ii+1][jj+1];  
	   }
     }
}

/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

void compute_max_warp_levels
(
                           /*************************************************/
 int   nx_orig,            /* in  : size in x-direction on finest grid      */
 int   ny_orig,            /* in  : size in y-direction on finest grid      */
 float n_eta,              /* in  : warping reduction factor                */
 int   *max_levels         /* out : maximum number of levels                */
                           /*************************************************/ 
)

/* compute maximum number of warping levels for given image size and warping */
/* reduction factor */

{

                           /*************************************************/
    int   i;               /* level counter                                 */
    float nx,ny;           /* reduced dimensions                            */
    float nx_old,ny_old;   /* reduced dimensions                            */
                           /*************************************************/

    nx_old=nx_orig;
    ny_old=ny_orig;

	/* Test different maximal number of warping levels i */
    for (i=1;;i++)
    { 
	  
	/* Compute corresponding resolution */
	nx=(int)ceil(nx_orig*pow(n_eta,i));
	ny=(int)ceil(ny_orig*pow(n_eta,i));	       

	/* if the resulting image size would be smaller than 4x4, we found the
	   maximal number of warp levels */
	if ((nx<4)||(ny<4)) break;

	nx_old=nx;
	ny_old=ny;
    } 
    
    if((nx==1)||(ny==1)) i--;
    *max_levels=i;
}



/* ------------------------------------------------------------------------- */

void update_nonlinearities_reg
(
 /*****************************************************/
 float **u,
 float **v,
 float **du,             /* in     : x-component of flow increment            */
 float **dv,             /* in     : y-component of flow increment            */
 float **psi_prime_s,    /* out    : nonlinearity data term                   */
 float lambda,        /* in     : diffusivity parameter                    */
 int   nx,               /* in     : size in x-direction                      */
 int   ny,               /* in     : size in y-direction                      */
 int   bx,               /* in     : boundary size in x-direction             */
 int   by,                /* in     : boundary size in y-direction             */
 float hx,
 float hy
/*****************************************************/
)
{
  int i,j;
  float help;
  float **ux, **uy, **vx, **vy;
  float **dux, **duy, **dvx, **dvy;
  malloc_multi(4,2,sizeof(float),nx,ny,bx,by,0,0,&ux,&uy,&vx,&vy);
  malloc_multi(4,2,sizeof(float),nx,ny,bx,by,0,0,&dux,&duy,&dvx,&dvy);
  
  mirror_bounds_2d(u, nx, ny, bx,by);
  mirror_bounds_2d(v, nx, ny, bx,by);
  mirror_bounds_2d(du, nx, ny, bx,by);
  mirror_bounds_2d(dv, nx, ny, bx,by);

  gradint_central( u, ux, uy, nx, ny, bx, by, hx, hy);
  gradint_central( v, vx, vy, nx, ny, bx, by, hx, hy);
  gradint_central( du, dux, duy, nx, ny, bx, by, hx, hy);
  gradint_central( dv, dvx, dvy, nx, ny, bx, by, hx, hy);

  
  
  
  for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {
      help = (ux[i][j]+dux[i][j])*(ux[i][j]+dux[i][j])
            +(vx[i][j]+dvx[i][j])*(vx[i][j]+dvx[i][j])
            +(uy[i][j]+duy[i][j])*(uy[i][j]+duy[i][j])
            +(vy[i][j]+dvy[i][j])*(vy[i][j]+dvy[i][j]);
      
      
      // Charbonnier function
      psi_prime_s[i][j] = 1.0/sqrt( help/(lambda * lambda) + 1.0);
    }
  
  
  free_multi(4,2,sizeof(float),nx,ny,bx,by,0,0,&ux,&uy,&vx,&vy);
  free_multi(4,2,sizeof(float),nx,ny,bx,by,0,0,&dux,&duy,&dvx,&dvy);
}



/* ------------------------------------------------------------------------- */

void update_nonlinearities
(
                        /*****************************************************/
float **J_11,           /* in     : entry 11 of the motion tensor            */
float **J_22,           /* in     : entry 22 of the motion tensor            */
float **J_33,           /* in     : entry 33 of the motion tensor            */
float **J_12,           /* in     : entry 12 of the motion tensor            */
float **J_13,           /* in     : entry 13 of the motion tensor            */
float **J_23,           /* in     : entry 23 of the motion tensor            */
float **du,             /* in     : x-component of flow increment            */
float **dv,             /* in     : y-component of flow increment            */
float **psi_prime_d,    /* out    : nonlinearity data term                   */
float epsilon_d,        /* in     : diffusivity parameter                    */
int   nx,               /* in     : size in x-direction                      */
int   ny,               /* in     : size in y-direction                      */
int   bx,               /* in     : boundary size in x-direction             */
int   by                /* in     : boundary size in y-direction             */
                        /*****************************************************/
)
{
int i,j;
float help;
for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {
		help =
					J_11[i][j]*du[i][j]*du[i][j]
			+       J_22[i][j]*dv[i][j]*dv[i][j]
			+ 2.0 * J_12[i][j]*du[i][j]*dv[i][j]
			+ 2.0 * J_13[i][j]*du[i][j]
			+ 2.0 * J_23[i][j]*dv[i][j]
			+       J_33[i][j];

			/* limited precision requires check for negative values */
			if (help<0) help=0;

      psi_prime_d[i][j] = 1.0/sqrt( help/(epsilon_d * epsilon_d) + 1.0);
	}
}
/* ------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------- */

void horn_schunck_warp_sor
(
 /*****************************************************/
 float **J_11,           /* in     : entry 11 of the motion tensor            */
 float **J_22,           /* in     : entry 22 of the motion tensor            */
 float **J_33,           /* in     : entry 33 of the motion tensor            */
 float **J_12,           /* in     : entry 12 of the motion tensor            */
 float **J_13,           /* in     : entry 13 of the motion tensor            */
 float **J_23,           /* in     : entry 23 of the motion tensor            */
 float **du,             /* in+out : x-component of flow increment            */
 float **dv,             /* in+out : y-component of flow increment            */
 float **u,              /* in     : x-component of flow field                */
 float **v,              /* in     : y-component of flow field                */
 float **psi_prime_d,    /* in     : nonlinearity data term                   */
 float **psi_prime_s,
 int   nx,               /* in     : size in x-direction                      */
 int   ny,               /* in     : size in y-direction                      */
 int   bx,               /* in     : boundary size in x-direction             */
 int   by,               /* in     : boundary size in y-direction             */
 float hx,               /* in     : grid spacing in x-direction              */
 float hy,               /* in     : grid spacing in y-direction              */
 float alpha,            /* in     : smoothness weight                        */
 float omega             /* in     : SOR overrelaxation parameter             */
/*****************************************************/
)

/*
 Computes one SOR iteration
 */

{
  /*****************************************************/
  int   i,j;              /* loop variables                                    */
  float hx_2,hy_2;        /* time saver variables                              */
  float xp,xm,yp,ym;      /* neighbourhood weights                             */
  float sum;              /* central weight                                    */
  /*****************************************************/
  
  
  /* define time saver variables */
  hx_2=alpha/(hx*hx);
  hy_2=alpha/(hy*hy);
  
  
  /* set boundaries zero */

  set_bounds_2d(du,nx,ny,bx,by,0.0);
  set_bounds_2d(dv,nx,ny,bx,by,0.0);
  
  /* mirror boundary	*/
  
  mirror_bounds_2d(u,nx,ny,bx,by);
  mirror_bounds_2d(v,nx,ny,bx,by);
  mirror_bounds_2d(psi_prime_s,nx,ny,bx,by);
  mirror_bounds_2d(psi_prime_d,nx,ny,bx,by);

  
  
  
  for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {
      /* compute weights */
      xp =  (i<nx+bx-1) * hx_2 * (psi_prime_s[i+1][j]+psi_prime_s[i][j])/2.0;
      xm =  (i>bx)      * hx_2 * (psi_prime_s[i-1][j]+psi_prime_s[i][j])/2.0;
      yp =  (j<ny+by-1) * hy_2 * (psi_prime_s[i][j+1]+psi_prime_s[i][j])/2.0;
      ym =  (j>by)      * hy_2 * (psi_prime_s[i][j-1]+psi_prime_s[i][j])/2.0;
      
      sum = xp + xm + yp + ym;
      
      /* perform SOR iteration */
      
      du[i][j]= (1.0-omega) * du[i][j] +
      omega *
      /* Gauss-Seidel step */
      (  - psi_prime_d[i][j]*J_13[i][j]
       - (psi_prime_d[i][j]*J_12[i][j] * dv[i  ][j  ]
          - xm  * (du[i-1][j  ] + u[i-1][j  ])
          - ym  * (du[i  ][j-1] + u[i  ][j-1])
          - yp  * (du[i  ][j+1] + u[i  ][j+1])
          - xp  * (du[i+1][j  ] + u[i+1][j  ])
          + sum * (               u[i  ][j  ])))
      /(psi_prime_d[i][j]*J_11[i][j]+sum);
      
      
      dv[i][j]= (1.0-omega) * dv[i][j] +
      omega *
      /* Gauss-Seidel step */
      (  - psi_prime_d[i][j]*J_23[i][j]
       - (psi_prime_d[i][j]*J_12[i][j] * du[i  ][j  ]
          - xm  * (dv[i-1][j  ] + v[i-1][j  ])
          - ym  * (dv[i  ][j-1] + v[i  ][j-1])
          - yp  * (dv[i  ][j+1] + v[i  ][j+1])
          - xp  * (dv[i+1][j  ] + v[i+1][j  ])
          + sum * (               v[i  ][j  ])))
      /(psi_prime_d[i][j]*J_22[i][j]+sum);
    }
  
}



/* ------------------------------------------------------------------------- */
void compute_motion_tensor
(
                        /*****************************************************/
float **f1,             /* in     : 1st image                                */
float **f2,             /* in     : 2nd image                                */
int   nx,               /* in     : size in x-direction                      */
int   ny,               /* in     : size in y-direction                      */
int   bx,               /* in     : boundary size in x-direction             */
int   by,               /* in     : boundary size in y-direction             */
float hx,               /* in     : grid spacing in x-direction              */
float hy,               /* in     : grid spacing in y-direction              */
float lambda,
float **J_11,           /* out    : entry 11 of the motion tensor            */
float **J_22,           /* out    : entry 22 of the motion tensor            */
float **J_33,           /* out    : entry 33 of the motion tensor            */
float **J_12,           /* out    : entry 12 of the motion tensor            */
float **J_13,           /* out    : entry 13 of the motion tensor            */
float **J_23            /* out    : entry 23 of the motion tensor            */
                        /*****************************************************/
)

/*
 Computes the motion tensor entries from a given image pair
*/

{
                        /*****************************************************/
int     i,j;            /* loop variables                                    */
float   **fx;           /* first order image derivatives                     */
float   **fy;           /* first order image derivatives                     */
float   **ft;           /* first order image derivatives                     */
float   **fxx;          /* second order image derivatives                    */
float   **fxy;          /* second order image derivatives                    */
float   **fyy;          /* second order image derivatives                    */
float   **fxt;          /* second order image derivatives                    */
float   **fyt;          /* second order image derivatives                    */
float   hx_1,hy_1;      /* time saver variables                              */
                        /*****************************************************/

float theta_x, theta_y, theta;// normalization factors for gradient-constant assumptions

/* allocate memory */

malloc_multi(8,2,sizeof(float),nx,ny,bx,by,0,0,&fx,&fy,&ft,&fxx,&fxy,&fyy,&fxt,&fyt);

/* define time saver variables */
hx_1=1.0/(2.0*hx);
hy_1=1.0/(2.0*hy);

/* mirror boundaries */
mirror_bounds_2d(f1,nx,ny,bx,by);
mirror_bounds_2d(f2,nx,ny,bx,by);

/* compute first oder derivatives */        
for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {
	fx[i][j] = 0.5*(f1[i+1][j]-f1[i-1][j]+f2[i+1][j]-f2[i-1][j])*hx_1;
	fy[i][j] = 0.5*(f1[i][j+1]-f1[i][j-1]+f2[i][j+1]-f2[i][j-1])*hy_1;
	ft[i][j] = (f2[i][j]-f1[i][j]);	
    }

/* mirror boundaries */
mirror_bounds_2d(fx,nx,ny,bx,by);
mirror_bounds_2d(fy,nx,ny,bx,by);
mirror_bounds_2d(ft,nx,ny,bx,by);

/* compute second order derivatives */
for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {
	fxx[i][j]=(fx[i+1][j  ]-fx[i-1][j  ])*hx_1;
	fxy[i][j]=(fy[i+1][j  ]-fy[i-1][j  ])*hx_1;
	fyy[i][j]=(fy[i  ][j+1]-fy[i  ][j-1])*hy_1;
	fxt[i][j]=(ft[i+1][j  ]-ft[i-1][j  ])*hx_1;
	fyt[i][j]=(ft[i  ][j+1]-ft[i  ][j-1])*hy_1;
    }


/* compute motion tensor entries entries */   
for(i=bx;i<nx+bx;i++)
    for(j=by;j<ny+by;j++)
    {  
 
	/* set up motion tensor */

	/* normalized Gradient constancy */

    theta_x = 1/(fxx[i][j] * fxx[i][j] + fxy[i][j] * fxy[i][j] + 0.1);
    theta_y = 1/(fxy[i][j] * fxy[i][j] + fyy[i][j] * fyy[i][j] + 0.1);
    theta = 1/(fx[i][j]*fx[i][j]+fy[i][j]*fy[i][j]+0.1);
//      theta_x = 1.0;
//      theta_y = 1.0;
//      theta = 1.0;
//	J_11[i][j] = theta_x*fxx[i][j] * fxx[i][j] + theta_y*fxy[i][j] * fxy[i][j];
//	J_22[i][j] = theta_x*fxy[i][j] * fxy[i][j] + theta_y*fyy[i][j] * fyy[i][j];
//	J_33[i][j] = theta_x*fxt[i][j] * fxt[i][j] + theta_y*fyt[i][j] * fyt[i][j];
//	J_12[i][j] = theta_x*fxx[i][j] * fxy[i][j] + theta_y*fxy[i][j] * fyy[i][j];
//	J_13[i][j] = theta_x*fxx[i][j] * fxt[i][j] + theta_y*fxy[i][j] * fyt[i][j];
//	J_23[i][j] = theta_x*fxy[i][j] * fxt[i][j] + theta_y*fyy[i][j] * fyt[i][j];
	  	
	
	/* Brightness constancy */
	/*
	J_11[i][j] = fx[i][j] * fx[i][j];
	J_22[i][j] = fy[i][j] * fy[i][j];
	J_33[i][j] = ft[i][j] * ft[i][j];
	J_12[i][j] = fx[i][j] * fy[i][j];
	J_13[i][j] = fx[i][j] * ft[i][j];
	J_23[i][j] = fy[i][j] * ft[i][j];
	*/

	/* Brightness constancy + 50* Gradient constancy */
//
      
	J_11[i][j] = lambda*(theta_x*fxx[i][j] * fxx[i][j] + theta_y*fxy[i][j]
                       * fxy[i][j]) + (1.0-lambda)*theta*fx[i][j] * fx[i][j];
	J_22[i][j] = lambda*(theta_x*fxy[i][j] * fxy[i][j] + theta_y*fyy[i][j]
                       * fyy[i][j]) + (1.0-lambda)*theta*fy[i][j] * fy[i][j];
	J_33[i][j] = lambda*(theta_x*fxt[i][j] * fxt[i][j] + theta_y*fyt[i][j]
                       * fyt[i][j]) + (1.0-lambda)*theta*ft[i][j] * ft[i][j];
	J_12[i][j] = lambda*(theta_x*fxx[i][j] * fxy[i][j] + theta_y*fxy[i][j]
                       * fyy[i][j]) + (1.0-lambda)*theta*fx[i][j] * fy[i][j];
	J_13[i][j] = lambda*(theta_x*fxx[i][j] * fxt[i][j] + theta_y*fxy[i][j]
                       * fyt[i][j]) + (1.0-lambda)*theta*fx[i][j] * ft[i][j];
	J_23[i][j] = lambda*(theta_x*fxy[i][j] * fxt[i][j] + theta_y*fyy[i][j]
                       * fyt[i][j]) + (1.0-lambda)*theta*fy[i][j] * ft[i][j];

    }




/* free memory */
free_multi(8,2,sizeof(float),nx,ny,bx,by,0,0,&fx,&fy,&ft,&fxx,&fxy,
           &fyy,&fxt,&fyt);
}

/* ------------------------------------------------------------------------- */

void HORN_SCHUNCK_WARP_LEVEL
(
                        /*****************************************************/
float **f1,             /* in     : 1st image                                */
float **f2,             /* in     : 2nd image                                */
float **du,             /* in+out : x-component of flow increment            */
float **dv,             /* in+out : y-component of flow increment            */
float **u,              /* in+out : x-component of flow field                */
float **v,              /* in+out : y-component of flow field                */
int   nx,               /* in     : size in x-direction                      */
int   ny,               /* in     : size in y-direction                      */
int   bx,               /* in     : boundary size in x-direction             */
int   by,               /* in     : boundary size in y-direction             */
float hx,               /* in     : grid spacing in x-direction              */
float hy,               /* in     : grid spacing in y-direction              */
float m_alpha,          /* in     : smoothness weight                        */
float epsilon_d,         /* in     : diffusivity param data term              */
float epsilon_s,       /* in     : diffusivity type data term               */
float w_bright_grad,
int   num_iter_inner,   /* in     : inner solver iterations                  */
int   num_iter_outer,   /* in     : outer nonlin update iterations           */
float n_omega           /* in     : SOR overrelaxation parameter             */

                        /*****************************************************/
)


{
                        /*****************************************************/
int    i,j;             /* loop variable                                     */
float  **J_11;          /* entry 11 of the motion tensor                     */
float  **J_22;          /* entry 22 of the motion tensor                     */
float  **J_33;          /* entry 33 of the motion tensor                     */
float  **J_12;          /* entry 12 of the motion tensor                     */
float  **J_13;          /* entry 13 of the motion tensor                     */
float  **J_23;          /* entry 23 of the motion tensor                     */
float  **psi_prime_d;   /* nonlinearitys of data term                        */
float  **psi_prime_s;
                        /*****************************************************/



/* ---- alloc memory ---- */
malloc_multi(8,2,sizeof(float),nx,ny,bx,by,0,0,&J_11,&J_22,&J_33,&J_12,&J_13,
             &J_23,&psi_prime_d, &psi_prime_s);

/* ---- initialise displacement field with zero ----  */
set_matrix_2d(du,nx+2*bx,ny+2*by,0,0,(float)0.0);
set_matrix_2d(dv,nx+2*bx,ny+2*by,0,0,(float)0.0);


/* ---- compute motion tensor ---- */
compute_motion_tensor(f1,f2,nx,ny,bx,by,hx,hy,w_bright_grad,
                      J_11,J_22,J_33,J_12,J_13,J_23);




/* ---- perform SOR iterations ---- */
for(j=0;j<num_iter_outer;j++)
{
	update_nonlinearities(J_11, J_22, J_33, J_12, J_13, J_23, du, dv,
                        psi_prime_d, epsilon_d, nx,
                        ny, bx, by);
  
  update_nonlinearities_reg(u,v,du,dv,psi_prime_s,epsilon_s,nx,ny,bx,by,
                            hx,hy);
	for(i=1;i<=num_iter_inner;i++)
	{
		horn_schunck_warp_sor(J_11, J_22, J_33, J_12, J_13, J_23,
							du, dv, u, v, psi_prime_d, psi_prime_s,nx, ny, bx, by, hx, hy,
							m_alpha, n_omega);
	}
  

}

/* ---- free memory ---- */
  free_multi(8,2,sizeof(float),nx,ny,bx,by,0,0,&J_11,&J_22,&J_33,&J_12,&J_13,
               &J_23,&psi_prime_d, &psi_prime_s);

}


/* ------------------------------------------------------------------------- */


void HORN_SCHUNCK_WARP
(
                        /*****************************************************/
float **f1_orig,        /* in     : 1st image (original resolution)          */
float **f2_orig,        /* in     : 2nd image (original resolution)          */
int   nx_orig,          /* in     : size in x-direction (original resolution)*/
int   ny_orig,          /* in     : size in y-direction (original resoluiton)*/
float **f1_res,         /* in+out : 1st image, resampled                     */
float **f2_res,         /* in+out : 2nd image, resampled                     */
float **f2_res_warp,    /* in+out : 2nd image, resampled  and warped         */
float **du,             /* in+out : x-component of flow increment            */
float **dv,             /* in+out : y-component of flow increment            */
float **u,              /* in+out : x-component of flow field                */
float **v,              /* in+out : y-component of flow field                */
float **tmp,            /* in+out : temporary aray for resampling            */
int   nx_fine,          /* in     : size in x-direction (current resolution) */
int   ny_fine,          /* in     : size in y-direction (current resolution) */
int   bx,               /* in     : boundary size in x-direction             */
int   by,               /* in     : boundary size in y-direction             */
float hx_fine,          /* in     : spacing in x-direction (current resol.)  */
float hy_fine,          /* in     : spacing in y-direction (current resol.)  */
float m_alpha,          /* in     : smoothness weight                        */
float epsilon_d,         /* in     : diffusivity param data term              */
float epsilon_s,
float w_bright_grad,
int   num_iter_inner,   /* in     : inner solver iterations                  */
int   num_iter_outer,   /* in     : outer nonlin update iterations           */
float n_omega,          /* in     : SOR overrelaxation parameter             */
float n_warp_eta,       /* in     : warping reduction factor between levels  */
int   max_rec_depth,    /* in     : maximum recursion depth                  */
int   rec_depth         /* in     : current recursion depth                  */

                        /*****************************************************/
)

/* implements warping for the Horn and Schunck method */

{

                             /************************************************/
int   nx_coarse,ny_coarse;   /* dimensions on previous coarser grid          */
float hx_coarse,hy_coarse;   /* grid sizes on previous coarser grid          */
                             /************************************************/
  

/* compute dimensions and grid sizes for previous coarser grid */	
nx_coarse=(int)ceil(nx_orig*pow(n_warp_eta,rec_depth+1));
ny_coarse=(int)ceil(ny_orig*pow(n_warp_eta,rec_depth+1));
hx_coarse=(float)nx_orig/(float)nx_coarse;
hy_coarse=(float)ny_orig/(float)ny_coarse;


/* start at coarsest level by recursively calling the routine */
if (rec_depth<max_rec_depth)
{
HORN_SCHUNCK_WARP(f1_orig, f2_orig, nx_orig, ny_orig,
		  f1_res, f2_res, f2_res_warp,
		  du, dv, u, v, tmp,
		  nx_coarse, ny_coarse, bx, by, hx_coarse, hy_coarse,	       
		  m_alpha, epsilon_d,epsilon_s, w_bright_grad,
		  num_iter_inner,num_iter_outer, n_omega, n_warp_eta, 
		  max_rec_depth,rec_depth+1);
}



/* ---- resample images ---------------------------------------------------- */

/* restrict original image pair to resolution of current level */ 
resample_2d(f1_orig,nx_orig,ny_orig,bx,by,f1_res,nx_fine,ny_fine,tmp);
resample_2d(f2_orig,nx_orig,ny_orig,bx,by,f2_res,nx_fine,ny_fine,tmp);


/* ---- get overall flow field from previous resolution level -------------- */

/* if on coarsest resolution */
if(rec_depth==max_rec_depth)
 {            
     /* set flow field zero */
     set_matrix_2d(u,nx_fine+2*bx,ny_fine+2*by,0,0,(float)0.0);
     set_matrix_2d(v,nx_fine+2*bx,ny_fine+2*by,0,0,(float)0.0); 
 }
 /* if not on coarsest resolution */
 else
 {
     /* interpolate solution from previous coarser level */    
     resample_2d(u,nx_coarse,ny_coarse,bx,by,u,nx_fine,ny_fine,tmp);
     resample_2d(v,nx_coarse,ny_coarse,bx,by,v,nx_fine,ny_fine,tmp);     
 }

/* ---- set up difference problem at current resolution -------------------- */
                
/* warp second image by overall flow field from previous coarser resolution */
backward_registration(f1_res,f2_res,f2_res_warp,u,v,
		              nx_fine,ny_fine,bx,by,hx_fine,hy_fine);



/* ---- solve difference problem at current resolution --------------------- */
     
/* solve difference problem at current resolution to obtain increment */
HORN_SCHUNCK_WARP_LEVEL(f1_res, f2_res_warp, du, dv, u, v,
                        nx_fine, ny_fine, bx, by, hx_fine, hy_fine,
                        m_alpha, epsilon_d, epsilon_s, w_bright_grad,
                        num_iter_inner, num_iter_outer, n_omega);

/* ---- compute overall flow field at current resolution ------------------- */
  
/* sum up flow increment */
add_matrix_2d(u,du,u,nx_fine,ny_fine,bx,by);
add_matrix_2d(v,dv,v,nx_fine,ny_fine,bx,by);
//  printf (" current depth: -- %d -- \n", rec_depth);
}



/* ------------------------------------------------------------------------- */


void HORN_SCHUNCK_MAIN
(
                         /*****************************************************/
float **f1,              /* in     : 1st image                                */
float **f2,              /* in     : 2nd image                                */
float **u,               /* out    : x-component of displacement field        */
float **v,               /* out    : y-component of displacement field        */
int   nx,                /* in     : size in x-direction                      */
int   ny,                /* in     : size in y-direction                      */
int   bx,                /* in     : boundary size in x-direction             */
int   by,                /* in     : boundary size in y-direction             */
float hx,                /* in     : grid spacing in x-direction              */
float hy,                /* in     : grid spacing in y-direction              */
float m_alpha,           /* in     : smoothness weight                        */
float epsilon_d,         /* in     : diffusivity param data term              */
float epsilon_s,
float w_bright_grad,
int   num_iter_inner,    /* in     : inner solver iterations                  */
int   num_iter_outer,    /* in     : outer nonlin update iterations           */
float n_warp_eta,        /* in     : warping reduction factor between levels  */
float n_omega,
int   n_warp_levels      /* in     : desired number of warping levels         */
                         /*****************************************************/
)

/* computes optic flow with Horn/Schunck + Warping */

{

                        /*****************************************************/
float **du;             /* x-component of flow increment                     */
float **dv;             /* y-component of flow increment                     */
float **f1_res;         /* 1st image, resampled                              */
float **f2_res;         /* 2nd image, resampled                              */
float **f2_res_warp;    /* 2nd image, resampled  and warped                  */
float **tmp;            /* temporary array for resampling                    */
int   max_rec_depth;    /* maximum recursion depth (warping level -1)        */
int   n_warp_max_levels;/* maximum possible number of warping levels         */
                        /*****************************************************/



/* compute maximal number of warping levels given the downsampling factor eta 
   and the image dimensions */
compute_max_warp_levels(nx,ny,n_warp_eta,&n_warp_max_levels);

/* limit number of desired warping levels by number of possible levels  */
/* we have to substract one, since: n warping levels -> n-1 recursions  */
max_rec_depth=minimum(n_warp_levels,n_warp_max_levels)-1;


/* ---- alloc memory ---- */

malloc_multi(6,2,sizeof(float),nx,ny,bx,by,0,0,&du,&dv,&f1_res,&f2_res,
             &f2_res_warp,&tmp);

/* call Horn/Schunck warping routine with desired number of levels */
HORN_SCHUNCK_WARP(f1, f2, nx, ny,
		  f1_res, f2_res, f2_res_warp,
		  du, dv, u, v, tmp,
		  nx, ny, bx, by, hx, hy,
		  m_alpha, epsilon_d, epsilon_s, w_bright_grad,
		  num_iter_inner, num_iter_outer, n_omega, n_warp_eta,
		  max_rec_depth,0);




/* ---- free memory ---- */
free_multi(6,2,sizeof(float),nx,ny,bx,by,0,0,&du,&dv,&f1_res,&f2_res,
           &f2_res_warp,&tmp);


}
/* ------------------------------------------------------------------------- */


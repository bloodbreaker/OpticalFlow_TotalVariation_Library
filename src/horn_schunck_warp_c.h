/*****************************************************************************/
/*                                                                           */
/*       Copyright 03/2012 by Dr. Andres Bruhn, Oliver Demetz and Yan Zhang  */
/*     Faculty of Mathematics and Computer Science, Saarland University,     */
/*                           Saarbruecken, Germany.                          */
/*                                                                           */
/*****************************************************************************/


#ifndef OF_HORN_SCHUNCK_INCLUDED
#define OF_HORN_SCHUNCK_INCLUDED


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
);

void AllocateMem2D( float ***img, int nx, int ny);
void DisallocateMem2D( float **img, int nx, int ny);





#endif

#ifndef NOISE_LIB_C
#define NOISE_LIB_C

#include <math.h>
#include <stdlib.h>
/*---------------------------------------------------------------------------*/
/*             copyright 2011 by Oliver Demetz                               */      
/*                               Mathematical Image Analysis Group           */
/*                               Saarland University                         */
/*---------------------------------------------------------------------------*/
void add_gauss_noise
(
	float **f,
	float **fn,
	float mean,
	float std,
	int nx,
	int ny,
	int bx,
	int by
)
{
	int i,j;
	float U,V,N,M;
	srand(1);
	for (i=bx; i<nx+bx; i++)
		for (j=by; j<ny+by; j++)
		{
			U = (float)rand() / RAND_MAX;
			V = (float)rand() / RAND_MAX;
			N= sqrtf(-2.0f*logf(U))*cosf(2.0f*M_PI*V);
			M= sqrtf(-2.0f*logf(U))*sinf(2.0f*M_PI*V);
			fn[i][j] = f[i][j] + N*std + mean;
		}
}

#endif

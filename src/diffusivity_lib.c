#ifndef DIFFUSIVITY_LIB
#define DIFFUSIVITY_LIB

/*---------------------------------------------------------------------------*/
/*             copyright 2011 by Oliver Demetz                               */      
/*                               Mathematical Image Analysis Group           */
/*                               Saarland University                         */
/*---------------------------------------------------------------------------*/

#include <math.h>
#include "funct_lib.c"

inline float pen_TV
(
	float s2,
	float epsilon
)
{
	return sqrtf(s2+epsilon*epsilon);
}

inline float diff_TV
(
	float s2,
	float epsilon
)
{
	return 1.0f / sqrtf(s2+epsilon*epsilon);
}

inline float diff_PM
(
	float s2,
	float lambda
)
{
	lambda *= lambda;
	return lambda / (lambda+s2);
}

inline float diff_Weickert
(
	float s2,
	float lambda
)
{
	return s2==0.0f ? 1.0f : 1.0f-expf(-3.31488/intpow(s2/(lambda*lambda),4));
}


inline float diff_Charbonnier
(
	float s2,
	float lambda
)
{
	return 1 / sqrtf(s2/(lambda*lambda)+1);
}

inline float diff_homo
(
	float s2,
	float lambda
)
{
	return 1.0f;
}

float diffusivity
(
	float s2,
	float lambda,
	int   diffusivity_type
)
{
	switch(diffusivity_type)
	{
		case 0: return diff_homo(s2,lambda);
		case 1: return diff_PM(s2,lambda);
		case 2: return diff_Weickert(s2,lambda);
		case 3: return diff_Charbonnier(s2,lambda);
		case 4: return diff_TV(s2,lambda);
		default: printf("invalid diffusivity_type specified\n");
	}
	return 0;
}

char diffusivity_names[5][20] = {"Homogeneous","PM","Weickert","Charbonnier","TV"};

#endif

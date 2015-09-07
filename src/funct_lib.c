/*****************************************************************************/
/*                                                                           */
/*         Copyright 08/2006 by Dr. Andres Bruhn and Oliver Demetz           */
/*     Faculty of Mathematics and Computer Science, Saarland University,     */
/*                           Saarbruecken, Germany.                          */
/*                                                                           */
/*****************************************************************************/


#ifndef FUNCT_LIB_INCLUDED
#define FUNCT_LIB_INCLUDED

#include <math.h>
#include <assert.h>

/*---------------------------------------------------------------------------*/
inline float maximumf(float a,float b)
{
	return a>b ? a : b;
}
/*---------------------------------------------------------------------------*/
inline float minimumf(float a,float b)
{
	return a<b ? a : b;
}
/*---------------------------------------------------------------------------*/
inline int maximum(int a,int b)
{
  return (a>b) ? a : b;
}
/*---------------------------------------------------------------------------*/
inline int minimum(int a,int b)
{
  return (a<b) ? a : b;
}
/*---------------------------------------------------------------------------*/
float intpow(float f, int exponent)
{
	if(exponent==0) return 1.0f;
	float result = f;
	int i;
	for (i=2;i<exponent;i*=2)
		result *= result;
	for ( i/=2;i<exponent;i++)
		result *=f;
	return result;
}
/*---------------------------------------------------------------------------*/
inline float square( float f)
{
	return f*f;
}
/*---------------------------------------------------------------------------*/
inline int isignf(float f)
{ 
	return f>=0.0f ? 1 : -1;
}
/*---------------------------------------------------------------------------*/
inline float sign(float f)
{
	if ( f>0.0f ) return 1.0f;
	else if ( f<0.0f ) return -1.0f;
	else return 0.0f;
}
/*---------------------------------------------------------------------------*/
float eucl_dist(float u1,float u2,float v1, float v2)
{
	return sqrtf(square(u1-v1)+square(u2-v2));
}
/*---------------------------------------------------------------------------*/
int count_occurences_in_string(const char* string, const char c)
{
	int i;
	int count=0;
	if(string!=0)
		for (i=0;i<strlen(string);i++)
		{
			if(string[i]==c) count++;
		}
	return count;
}
/*---------------------------------------------------------------------------*/
inline void clamp (int min, int *val, int max)
{
	assert(min<=max);
	
	if (*val<min) *val=min;
	if (*val>max) *val=max;
}
/*---------------------------------------------------------------------------*/
inline void clampf (float min, float *val, float max)
{
	assert(min<=max);
	
	if (*val<min) *val=min;
	if (*val>max) *val=max;
}
/*---------------------------------------------------------------------------*/
inline int clamped (int min, int val, int max)
{
	assert(min<=max);
	return minimum(maximum(min,val),max);
}

inline float clampedf (float min, float val, float max)
{
	assert(min<=max);
	return minimumf(maximumf(min,val),max);
}
/*---------------------------------------------------------------------------*/
int  in_range (int val, int mn, int mx)
{
	return (val >= mn && val < mx);
}
/*---------------------------------------------------------------------------*/
int  in_range_f (float val, float mn, float mx)
{
	return (val >= mn && val < mx);
}
/*---------------------------------------------------------------------------*/
int  byte_range(int a)
/* restricts number to unsigned char range */
{
    return clamped(0,a,255);
}
/*---------------------------------------------------------------------------*/
unsigned char float_to_uchar (float f)
{
	return (unsigned char) clamped(0,(int)(f+0.5f),255);
}
/*---------------------------------------------------------------------------*/
#endif

#ifndef FRONTEND_TOOLS
#define FRONTEND_TOOLS

/*---------------------------------------------------------------------------*/
/*             copyright 2011 by Oliver Demetz                               */      
/*                               Mathematical Image Analysis Group           */
/*                               Saarland University                         */
/*---------------------------------------------------------------------------*/


#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*---------------------------------------------------------------------------*/

void change_vali(	int key,
					unsigned char active_param,
					const unsigned char assigned_param,
					int *val)
{
	if (active_param==assigned_param)
	{
		if(key==GLUT_KEY_UP)
			*val = *val + 1;
		else if (key==GLUT_KEY_DOWN)
			*val = *val - 1;
	}
}

/*---------------------------------------------------------------------------*/

void change_vali_cyclic(	int key,
							unsigned char active_param,
							const unsigned char assigned_param,
							int *val, 
							int mod)
{
	if (active_param==assigned_param)
	{
		if(key==GLUT_KEY_UP)
			*val = (*val + 1)%mod;
		else if (key==GLUT_KEY_DOWN)
			*val = (*val - 1)%mod;
	}
}

/*---------------------------------------------------------------------------*/

void change_valf(	int key,
					unsigned char active_param,
					const unsigned char assigned_param,
					float* val)
{
	if (active_param==assigned_param)
	{
		if(key==GLUT_KEY_UP)
		{
			if(*val==0.0f)
				*val = 0.01f;
			else
				*val = (*val) * 1.1f;
		}
		else if (key==GLUT_KEY_DOWN)
		{
			if(*val>0.001f)
				*val = *val / 1.1f;
			else
				*val = 0.0f;
		}
	}
}

/*---------------------------------------------------------------------------*/

void showParamLinei(	char * str,
						unsigned char active_param,
						const unsigned char assigned_key,
						const char* param_long_name,
						int *val)
{
	char temp[100];
	sprintf(temp,"%c %s %c",active_param==assigned_key?'(':' ',
									param_long_name,
									active_param==assigned_key?')':' ');
	sprintf(str+strlen(str)," (%c) %40s   %12d\n",	assigned_key,temp,*val);
}

/*---------------------------------------------------------------------------*/

void showParamLinef(	char * str,
						unsigned char active_param,
						const unsigned char assigned_key,
						const char* param_long_name,
						float *val)
{
	char temp[100];
	sprintf(temp,"%c %s %c",active_param==assigned_key?'(':' ',
									param_long_name,
									active_param==assigned_key?')':' ');
	sprintf(str+strlen(str)," (%c) %40s   %12.4f\n",	assigned_key,temp,*val);
}

/*---------------------------------------------------------------------------*/

void showParamLineStr(	char * str,
						unsigned char active_param,
						const unsigned char assigned_key,
						const char* param_long_name,
						char *val)
{
	char temp[100];
	sprintf(temp,"%c %s %c",active_param==assigned_key?'(':' ',
									param_long_name,
									active_param==assigned_key?')':' ');
	sprintf(str+strlen(str)," (%c) %40s   %12s\n",	assigned_key,temp,val);
}

/*---------------------------------------------------------------------------*/

int __index=0;

char names[5][20] = {	"lastcompute0.pgm",
						"lastcompute1.pgm",
						"lastcompute2.pgm",
						"lastcompute3.pgm",
						"lastcompute4.pgm"};
char* get_lastname()
{
	return names[(__index++)%5];
}

void printTitleLine(	char * str,
						const char* title)
{
int i;
int len = strlen(title);
int tgtlen=60;
for(i=0;i<(tgtlen-len-2)/2;i++)
	sprintf(str+strlen(str),"-");
sprintf(str+strlen(str)," %s ", title);
for(i=0;i<(tgtlen-len-2)/2;i++)
	sprintf(str+strlen(str),"-");
sprintf(str+strlen(str),"\n");
}

void printFooter(	char * str)
{
int i;
for(i=0;i<60;i++)
	sprintf(str+strlen(str),"-");
sprintf(str+strlen(str),"\n");
}

#endif
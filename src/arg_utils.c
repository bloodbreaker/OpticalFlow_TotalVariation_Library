#ifndef ARG_UTILS
#define ARG_UTILS

#include <stdio.h>

void parse_arg_float ( int argc, char** argv, const char* argname, float* variable,float stdval, char* reportstring)
{
	int i;
	int found=0;
	if (argc%2!=1)
	{
		printf("Invalid number of arguments! Exiting.\n");
		exit(1);
	}
	
	for (i=1;i<argc-1;i++)
	{
		if(strcmp(argv[i],argname)==0)
		{
			*variable = atof(argv[i+1]);
			found=1;
			break;
		}
	}
	if(found==0)
	{
		*variable = stdval;
	}
	sprintf(reportstring+strlen(reportstring),"f %s = %f %s\n", argname, *variable,(found==0 ? "" :"*"));
}

void parse_arg_int ( int argc, char** argv, const char* argname, int* variable, int stdval, char* reportstring)
{
	int i;
	int found=0;
	if (argc%2!=1)
	{
		printf("Invalid number of arguments! Exiting.\n");
		exit(1);
	}
	
	for (i=1;i<argc-1;i++)
	{
		if(strcmp(argv[i],argname)==0)
		{
			*variable = atoi(argv[i+1]);
			found=1;
			break;
		}
	}
	if(found==0)
	{
		*variable = stdval;
	}
	sprintf(reportstring+strlen(reportstring),"i %s = %d %s\n", argname, *variable,(found==0 ? "" :"*"));
}

void parse_arg_string ( int argc, char** argv, const char* argname, char* variable,const char* stdval, char* reportstring)
{
	int i;
	int found=0;
	if (argc%2!=1)
	{
		printf("Invalid number of arguments! Exiting.\n");
		exit(1);
	}
	
	for (i=1;i<argc-1;i++)
	{
		if(strcmp(argv[i],argname)==0)
		{
			sprintf(variable,"%s",argv[i+1]);
			found=1;
		}
	}
	if(found==0)
	{
		sprintf(variable,"%s",stdval);
	}
	sprintf(reportstring+strlen(reportstring),"s %s = %s %s\n", argname, variable,(found==0 ? "" :"*"));
}

#endif
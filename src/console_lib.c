#ifndef CONSOLE_LIB
#define CONSOLE_LIB
/*****************************************************************************/
/* --- console_lib --------------------------------------------------------- */
/* Output console stuff                                                      */
/*                                                                           */
/* (C) 2008 - 2010 Pascal Gwosdek                                            */
/*                 Mathematical Image Analysis Group,                        */
/*                 Saarland University, Germany                              */
/*                                                                           */
/* Version 1.01 (2010-09-15)                                                 */
/*****************************************************************************/

#ifndef NDEBUG
#define DEBUG_OUT(...) console_debug(__FILE__,__LINE__,__VA_ARGS__)
#else
#define DEBUG_OUT(...)
#endif

#include <string.h>
#include <stdlib.h>
#include <stdarg.h>

FILE * shadowfile = 0;

void init_shadow(const char *formatstring, ...)
{
	va_list args;
	char fname[255];
	va_start(args, formatstring);
	vsprintf(fname, formatstring, args);
	if(shadowfile) fclose(shadowfile);
	shadowfile = fopen(fname,"w");
	va_end(args);
}

/*****************************************************************************/
/* Outputs a message.                                                        */
/*****************************************************************************/
void console_out(const char *formatstring, ...)
{
  va_list args;
  va_start(args, formatstring);
  vfprintf(stdout, formatstring, args);
  if(shadowfile)
  {
	  va_end(args);
	  va_start(args, formatstring);
	  vfprintf(shadowfile, formatstring, args);
  }
  
  int pos = strlen(formatstring);

  if (formatstring[pos-1] != '\n')
  {
    fprintf(stdout, "\n");
	if(shadowfile)
		fprintf(shadowfile, "\n");
  }
  if(shadowfile) fflush(shadowfile);
  va_end(args);
}

/*****************************************************************************/
/* Outputs a message with a green INFO label.                                */
/*****************************************************************************/
void console_info(const char *formatstring, ...)
{
  va_list args;
  va_start(args, formatstring);
  fprintf(stdout, "\033[32;1mINFO:\033[0m    ");
  vfprintf(stdout, formatstring, args);

  int pos = strlen(formatstring);

  if (formatstring[pos-1] != '\n')
    fprintf(stdout, "\n");

  va_end(args);
}

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Outputs a message with a red ERROR label.                                 */
/*****************************************************************************/
void console_error(const char *formatstring, ...)
{
  va_list args;
  va_start(args, formatstring);
  fprintf(stderr, "\033[31;1mERROR:\033[0m   ");
  vfprintf(stderr, formatstring, args);

  int pos = strlen(formatstring);

  if (formatstring[pos-1] != '\n')
    fprintf(stderr, "\n");

  va_end(args);
}

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Outputs a message with a yellow WARNING label.                            */
/*****************************************************************************/
void console_warning(const char *formatstring, ...)
{
  va_list args;
  va_start(args, formatstring);
  fprintf(stderr, "\033[33;1mWARNING:\033[0m ");
  vfprintf(stderr, formatstring, args);

  int pos = strlen(formatstring);

  if (formatstring[pos-1] != '\n')
    fprintf(stderr, "\n");

  va_end(args);
}

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Outputs a message with a blue DEBUG label.                                */
/*****************************************************************************/
void console_debug(const char *file, int line, const char *formatstring, ...)
{
  va_list args;
  va_start(args, formatstring);
  fprintf(stderr, "\033[34;1mDEBUG:\033[0m   ");
  fprintf(stderr, "(%s:%d): ", file, line);
  vfprintf(stderr, formatstring, args);

  int pos = strlen(formatstring);

  if (formatstring[pos-1] != '\n')
    fprintf(stderr, "\n");

  fflush(stderr);

  va_end(args);
}

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Outputs a message with a enough space to the left, to continue a message  */
/*****************************************************************************/
void console_text(const char *formatstring, ...)
{
  va_list args;
  va_start(args, formatstring);
  fprintf(stdout, "         ");
  vfprintf(stdout, formatstring, args);

  int pos = strlen(formatstring);

  if (formatstring[pos-1] != '\n')
    fprintf(stdout, "\n");

  va_end(args);
}

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Outputs a message with a enough space to the left, to continue a message  */
/* different to _info                                                        */
/*****************************************************************************/
void console_text_noinfo(const char *formatstring, ...)
{
  va_list args;
  va_start(args, formatstring);
  fprintf(stderr, "         ");
  vfprintf(stderr, formatstring, args);

  int pos = strlen(formatstring);

  if (formatstring[pos-1] != '\n')
    fprintf(stderr, "\n");

  va_end(args);
}

#endif

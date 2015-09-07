/*****************************************************************************/
/* --- io_lib -------------------------------------------------------------- */
/* Input and output images                                                   */
/*                                                                           */
/* (C) 2008 - 2011 Pascal Gwosdek and Oliver Demetz                          */
/*                 Mathematical Image Analysis Group,                        */
/*                 Saarland University, Germany                              */
/*                                                                           */
/* Version 1.04 (2011-06-06)                                                 */
/*****************************************************************************/
#ifndef IO_LIB_C
#define IO_LIB_C

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include "malloc_lib.c"
#include "console_lib.c"
#include "funct_lib.c"
/*****************************************************************************/
/* Image types                                                               */
/*****************************************************************************/

typedef enum
{
  IMG_TYPE_INVALID = 0,
  IMG_TYPE_PGM     = 2,
  IMG_TYPE_PPM     = 3
} image_t;

#include <stdarg.h>

/*****************************************************************************/
/* Read in a one-byte or two-byte PGM file.                                  */
/* RETURNS a pointer to the float-based image contents with bx,by wide       */
/*         boundaries (uninitialised), or NULL on an error.                  */
/*                                                                           */
/* The memory is allocated with sufficient padding, such that width and      */
/* height can be addressed up to the next multiple of padding_x and          */
/* padding_y, respectively. If bx and by are specified additionally, they    */
/* will additionally be added to the padded result.                          */
/*                                                                           */
/* The allocated memory needs to be freed by the user program afterwards,    */
/* by calling free_pgm_image (padding is irrelevant for 'free').             */
/* On success, the width and height of the image are in width and height.    */
/*****************************************************************************/
float **read_pgm_image_padded
(
  const char    *filename,      /* Name of file to be read                   */
  int           *width,         /* Returned image width                      */
  int           *height,        /* Returned image height                     */
  int           bx,             /* Boundary in X                             */
  int           by,             /* Boundary in Y                             */
  int           *max_grey_value,/* Maximal grey value                        */
  int           padding_x,      /* Padding to that many elements in x dir.   */
  int           padding_y,      /* Padding to that many elements in y dir.   */
  int           *memwidth,      /* Return: Resulting #elements in memory (x) */
  int           *memheight      /* Return: Resulting #elements in memory (y) */
)
{
  FILE         *image_file;
  char         *line = NULL;
  float        **image;
  size_t       length = 0;

  if (!(image_file = fopen(filename, "r")))
  {
    *width    = 0; *height    = 0;
    *memwidth = 0; *memheight = 0;
    console_error("Could not open file %s!\n", filename);
    return NULL;
  }

  /* Check for the magic bytes                                               */
  if (getline(&line, &length, image_file) >= 2)
  {
    if ((line[0] != 'P') || (line[1] != '5'))
    {
      *width    = 0; *height    = 0;
      *memwidth = 0; *memheight = 0;
      console_error("File %s if not of type PGM!\n", filename);
      fclose(image_file);
      return NULL;
    }
  }
  else
  {
    *width    = 0; *height    = 0;
    *memwidth = 0; *memheight = 0;
    console_error("Unexpected end of file %s!\n", filename);
    fclose(image_file);
    return NULL;
  }

  /* Read in the image dimensions                                            */
  *width    = 0; *height    = 0;
  *memwidth = 0; *memheight = 0;

  while (getline(&line, &length, image_file) > 0)
  {
    if((line[0] != '#') && (line[0] != '\n'))
    {
      sscanf(line, "%d %d", width, height);
      break;
    }
  }

  if ((*width == 0) || (*height == 0))
  {
    console_error("Unexpected end of file %s!\n", filename);
    fclose(image_file);
    return NULL;
  }

  /* Account for padding                                                             */
  if ((padding_x != 1) || (padding_y != 1))
  {
    if (!*width || !*height)
    {
      console_error("Please specify meaningful paddings (default: 1)!");
      *width    = 0; *height    = 0;
      *memwidth = 0; *memheight = 0;
      return NULL;
    }
    
    *memwidth  = ((*width  - 1) / padding_x + 1) * padding_x;
    *memheight = ((*height - 1) / padding_y + 1) * padding_y;
  }
  else
  {
    *memwidth  = *width;
    *memheight = *height;
  }

  /* Read in the image *max_grey_value                                               */
  *max_grey_value = 0;
  while (getline(&line, &length, image_file) > 0)
  {
    if((line[0] != '#') && (line[0] != '\n'))
    {
      sscanf(line, "%d", &*max_grey_value);
      break;
    }
  }

  if (*max_grey_value == 0)
  {
    *width    = 0; *height    = 0;
    *memwidth = 0; *memheight = 0;
    console_error("Unexpected end of file %s!\n", filename);
    fclose(image_file);
    return NULL;
  }

  if (line != NULL)
  {
    free(line); line = NULL;
  }

  /* Allocate the image array                                                */
  calloc_multi(1,         2,          sizeof(float),
               *memwidth, *memheight,
               bx,        by,
               0,         0,
               &image);

  
  /* Read in the raw input                                                   */
  if (*max_grey_value < 256)
  {
    int number = 0;
    int i, j;

    for (j = by; j < *height + by; ++j)
      for (i = bx; i < *width + bx; ++i)
      {
        if ((number = fgetc(image_file)) != EOF)
        {
          image[i][j] = (float)number;
        }
        else
        {
          *width    = 0; *height    = 0;
          *memwidth = 0; *memheight = 0;
          free_multi(1,         2,          sizeof(float),
                     *memwidth, *memheight,
                     bx,        by,
                     0,         0,
                     &image);
          console_error("Unexpected end of file %s!\n", filename);
          fclose(image_file);
          return NULL;
        }
      }
  }
  else if (*max_grey_value < 65536)
  {
    unsigned int number = 0;
    int          i, j;

    for (j = by; j < *height + by; ++j)
      for (i = bx; i < *width + bx; ++i)
      {
        if (fscanf(image_file, "%2c",
                   (char*)&number + sizeof(unsigned int) - 2) != EOF)
        {
          image[i][j] = (float)number;
        }
        else
        {
          *width    = 0; *height    = 0;
          *memwidth = 0; *memheight = 0;
          free_multi(1,         2,          sizeof(float),
                     *memwidth, *memheight,
                     bx,        by,
                     0,         0,
                     &image);
          console_error("Unexpected end of file %s!\n", filename);
          fclose(image_file);
          return NULL;
        }
      }
  }
  else
  {
    *width    = 0; *height    = 0;
    *memwidth = 0; *memheight = 0;
    console_error("Unexpected max grey value in file %s!\n", filename);
    fclose(image_file);
    return NULL;
  }

  fclose(image_file);
  return image;
}

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Read in a one-byte or two-byte PGM file.                                  */
/* RETURNS a pointer to the float-based image contents with bx,by wide       */
/*         boundaries (uninitialised), or NULL on an error.                  */
/*                                                                           */
/* The allocated memory needs to be freed by the user program afterwards,    */
/* either manually or by calling free_pgm_image.                             */
/* On success, the width and height of the image are in width and height.    */
/*****************************************************************************/
float **read_pgm_image
(
  const char    *filename,      /* Name of file to be read                   */
  int           *width,         /* Returned image width                      */
  int           *height,        /* Returned image height                     */
  int           bx,             /* Boundary in X                             */
  int           by,             /* Boundary in Y                             */
  int           *max_grey_value /* Maximal grey value                        */
)
{
  int    memwidth, memheight;
  return read_pgm_image_padded(filename, width, height, bx, by, max_grey_value,
                               1, 1, &memwidth, &memheight);
}

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Oliver Demetz: June 2011                                                  */
/* Internal function that processes the va_list                              */
/*****************************************************************************/
void write_pgm_image_comment_started
(
  const char    *filename,      /* Name of file to be written                */
  int           width,          /* Width of the image                        */
  int           height,         /* Height of the image                       */
  int           bx,             /* Boundary in X                             */
  int           by,             /* Boundary in Y                             */
  float         **image,        /* Input image                               */
  int           bytes,          /* Byte format to use [1,2]                  */
  const char    *comment,       /* Comment lines                             */
  va_list       args            /* variable length argument list             */
)
{
  FILE         *image_file;
  va_list      args_copy;
  int          written;
  char         *filename_replaced;
  
  va_copy(args_copy,args);
  
  if (image == 0)
  {
    console_error("Could not write image %s - NULL pointer! Maybe there\n",
                  filename);
    console_text ("was a problem loading the image already?\n");
    return;
  }

  if ((bytes < 1) || (bytes > 2))
  {
    console_error("Unexpected format for file %s!\n", filename);
    return;
  }

  filename_replaced = (char*)calloc(1024,sizeof(char));
  written = vsnprintf (filename_replaced,1024,filename, args);
  if (!(image_file = fopen(filename_replaced, "w")))
  {
    console_error("Could not open file %s!\n", filename_replaced);
    return;
  }

  /* Write magic bytes                                                       */
  fprintf(image_file, "P5\n");

  /* Write (optional) comment lines                                          */
  if (comment)
  {
	int pos     = 0;
    int newline = 1;
    char *fullstring;
    char *comment_replaced;
	char *format;
	
	format = (char*)calloc(strlen(filename)+strlen(comment)+1,sizeof(char));
	strcpy(format,filename);
	strcat(format,comment);
	fullstring = (char*)calloc(2048, sizeof(char));
    vsnprintf (fullstring, 2048, format, args_copy);
	comment_replaced=fullstring+written;
	
    while (comment_replaced[pos] != '\0')
    {
      if (newline && (comment_replaced[pos] != '#'))
      {
        fprintf(image_file, "#%c", comment_replaced[pos]);
        newline = 0;
      }
      else
        fprintf(image_file, "%c", comment_replaced[pos]);

      if (comment_replaced[pos] == '\n')
        newline = 1;

      ++pos;
    }

    if (!newline)
      fprintf(image_file, "\n");
	
	free(fullstring);free(format);
  }

  /* Write dimensions                                                        */
  fprintf(image_file, "%d %d\n", width, height);

  /* Write depth                                                             */
  if (bytes == 1)
    fprintf(image_file, "255\n");
  else if (bytes == 2)
    fprintf(image_file, "65535\n");

  /* Read in the raw input                                                   */
  if (bytes == 1)
  {
    int i, j;

    for (j = by; j < height + by; ++j)
      for (i = bx; i < width + bx; ++i)
        fprintf(image_file, "%c",float_to_uchar(image[i][j]));
  }
  else if (bytes == 2)
  {
    int i, j;

    for (j = by; j < height + by; ++j)
      for (i = bx; i < width + bx; ++i)
        fprintf(image_file, "%2c",float_to_uchar(image[i][j]));
  }

  fclose(image_file);
  va_end (args_copy);
  printf("written %s\n",filename_replaced);
  free(filename_replaced);
}

/*****************************************************************************/
/* Write a one-byte or two-byte PGM file.                                    */
/*                                                                           */
/* Comment takes several comment lines, separated by '\n'. There is no need  */
/* to include markers such as '#', since they are automatically inserted     */
/* (where necessary)                                                         */
/* O.Demetz, June 2011: Extended functionality of printf-like features for   */
/* the filename and comment. Just give filename and comment with %s and %f   */
/* and such patterns, and append the corresponding variables like for printf.*/
/*****************************************************************************/
void write_pgm_image_comment
(
  const char    *filename,      /* Name of file to be written                */
  int           width,          /* Width of the image                        */
  int           height,         /* Height of the image                       */
  int           bx,             /* Boundary in X                             */
  int           by,             /* Boundary in Y                             */
  float         **image,        /* Input image                               */
  int           bytes,          /* Byte format to use [1,2]                  */
  const char    *comment,       /* Comment lines                             */
  ...
)
{
	va_list args;
	va_start(args,comment);
	write_pgm_image_comment_started(filename,width,height,bx,by,image,bytes,
									comment,args);
	va_end(args);
}



/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Write a one-byte or two-byte PGM file.                                    */
/*                                                                           */
/* Comment takes several comment lines, separated by '\n'. There is no need  */
/* to include markers such as '#', since they are automatically inserted     */
/* (where necessary)                                                         */
/* O.Demetz, June 2011: Extended functionality of printf-like features for   */
/* the filename and comment. Just give filename and comment with %s and %f   */
/* and such patterns, and append the corresponding variables like for printf.*/
/*****************************************************************************/
void write_pgm_image
(
  const char    *filename,      /* Name of file to be written                */
  int           width,          /* Width of the image                        */
  int           height,         /* Height of the image                       */
  int           bx,             /* Boundary in X                             */
  int           by,             /* Boundary in Y                             */
  float         **image,        /* Input image                               */
  int           bytes,          /* Byte format to use [1,2]                  */
  ...
)
{
  va_list args;
  va_start(args,bytes);
  write_pgm_image_comment_started(filename, width, height, bx, by, 
								  image, bytes, NULL, args);
  va_end(args);
}

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Read in a one-byte or two-byte PPM file.                                  */
/* RETURNS a pointer to the float-based image contents with bx,by wide       */
/*         boundaries (uninitialised), or NULL on an error.                  */
/*                                                                           */
/* The memory is allocated with sufficient padding, such that width and      */
/* height can be addressed up to the next multiple of padding_x and          */
/* padding_y, respectively. If bx and by are specified additionally, they    */
/* will additionally be added to the padded result.                          */
/*                                                                           */
/* The allocated memory needs to be freed by the user program afterwards,    */
/* either manually or by calling free_ppm_image.                             */
/* On success, the width and height of the image are in width and height.    */
/*****************************************************************************/
float ***read_ppm_image_padded
(
  const char    *filename,      /* Name of file to be read                   */
  int           *width,         /* Returned image width                      */
  int           *height,        /* Returned image height                     */
  int           bx,             /* Boundary in X                             */
  int           by,             /* Boundary in Y                             */
  int           *max_grey_value,/* Maximal grey value                        */
  int           padding_x,      /* padding to that many elements in x dir.   */
  int           padding_y,      /* padding to that many elements in y dir.   */
  int           *memwidth,      /* Return: Resulting #elements in memory (x) */
  int           *memheight      /* Return: Resulting #elements in memory (y) */
)
{
  FILE         *image_file;
  char         *line = NULL;
  float        ***image;
  size_t       length = 0;

  if (!(image_file = fopen(filename, "r")))
  {
    *width    = 0; *height    = 0;
    *memwidth = 0; *memheight = 0;
    console_error("Could not open file %s!\n", filename);
    return NULL;
  }

  /* Check for the magic bytes                                               */
  if (getline(&line, &length, image_file) >= 2)
  {
    if ((line[0] != 'P') || (line[1] != '6'))
    {
      *width    = 0; *height    = 0;
      *memwidth = 0; *memheight = 0;
      console_error("File %s if not of type PPM!\n", filename);
      fclose(image_file);
      return NULL;
    }
  }
  else
  {
    *width    = 0; *height    = 0;
    *memwidth = 0; *memheight = 0;
    console_error("Unexpected end of file %s!\n", filename);
    fclose(image_file);
    return NULL;
  }

  /* Read in the image dimensions                                            */
  *width    = 0; *height    = 0;
  *memwidth = 0; *memheight = 0;

  while (getline(&line, &length, image_file) > 0)
  {
    if((line[0] != '#') && (line[0] != '\n'))
    {
      sscanf(line, "%d %d", width, height);
      break;
    }
  }

  if ((*width == 0) || (*height == 0))
  {
    console_error("Unexpected end of file %s!\n", filename);
    fclose(image_file);
    return NULL;
  }

  /* Account for padding                                                             */
  if ((padding_x != 1) || (padding_y != 1))
  {
    if (!*width || !*height)
    {
      console_error("Please specify meaningful paddings (default: 1)!");
      *width    = 0; *height    = 0;
      *memwidth = 0; *memheight = 0;
      return NULL;
    }
    
    *memwidth  = ((*width  - 1) / padding_x + 1) * padding_x;
    *memheight = ((*height - 1) / padding_y + 1) * padding_y;
  }
  else
  {
    *memwidth  = *width;
    *memheight = *height;
  }

  /* Read in the image *max_grey_value                                               */
  *max_grey_value = 0;
  while (getline(&line, &length, image_file) > 0)
  {
    if((line[0] != '#') && (line[0] != '\n'))
    {
      sscanf(line, "%d", &*max_grey_value);
      break;
    }
  }

  if (*max_grey_value == 0)
  {
    *width    = 0; *height    = 0;
    *memwidth = 0; *memheight = 0;
    console_error("Unexpected end of file %s!\n", filename);
    fclose(image_file);
    return NULL;
  }

  if (line != NULL)
  {
    free(line); line = NULL;
  }

  /* Allocate the image array                                                */
  malloc_multi(1, 3,         sizeof(float),
               3, *memwidth, *memheight,
               0, bx,        by,
               0, 0,         0,
               &image);
  
  /* Read in the raw input                                                   */
  if (*max_grey_value < 256)
  {
    int number = 0;
    int i, j, c;

    for (j = by; j < *height + by; ++j)
      for (i = bx; i < *width + bx; ++i)
        for (c = 0; c < 3; ++c)
        {
          if ((number = fgetc(image_file)) != EOF)
          {
            image[c][i][j] = (float)number;
          }
          else
          {
            *width    = 0; *height    = 0;
            *memwidth = 0; *memheight = 0;
            free_multi(1, 3,         sizeof(float),
                       3, *memwidth, *memheight,
                       0, bx,        by,
                       0, 0,         0,
                       &image);
            console_error("Unexpected end of file %s!\n", filename);
            fclose(image_file);
            return NULL;
          }
        }
  }
  else if (*max_grey_value < 65536)
  {
    unsigned int number = 0;
    int          i, j, c;

    for (j = by; j < *height + by; ++j)
      for (i = bx; i < *width + bx; ++i)
        for (c = 0; c < 3; ++c)
        {
          if (fscanf(image_file, "%2c",
                     (char*)&number + sizeof(unsigned int) - 2) != EOF)
          {
            image[c][i][j] = (float)number;
          }
          else
          {
            *width    = 0; *height    = 0;
            *memwidth = 0; *memheight = 0;
            free_multi(1, 3,          sizeof(float),
                       3, *memwidth, *memheight,
                       0, bx,        by,
                       0, 0,         0,
                       &image);
            console_error("Unexpected end of file %s!\n", filename);
            fclose(image_file);
            return NULL;
          }
        }
  }
  else
  {
    *width    = 0; *height    = 0;
    *memwidth = 0; *memheight = 0;
    console_error("Unexpected max grey value in file %s!\n", filename);
    fclose(image_file);
    return NULL;
  }

  fclose(image_file);
  return image;
}

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Read in a one-byte or two-byte PPM file.                                  */
/* RETURNS a pointer to the float-based image contents with bx,by wide       */
/*         boundaries (uninitialised), or NULL on an error.                  */
/*                                                                           */
/* The allocated memory needs to be freed by the user program afterwards,    */
/* either manually or by calling free_ppm_image.                             */
/* On success, the width and height of the image are in width and height.    */
/*****************************************************************************/
float ***read_ppm_image
(
  const char    *filename,      /* Name of file to be read                   */
  int           *width,         /* Returned image width                      */
  int           *height,        /* Returned image height                     */
  int           bx,             /* Boundary in X                             */
  int           by,             /* Boundary in Y                             */
  int           *max_grey_value /* Maximal grey value                        */
)
{
  int    memwidth, memheight;
  return read_ppm_image_padded(filename, width, height, bx, by, max_grey_value,
                               1, 1, &memwidth, &memheight);
};

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Oliver Demetz: June 2011                                                  */
/* Internal function that processes the va_list                              */
/*****************************************************************************/
void write_ppm_image_comment_started
(
  const char    *filename,      /* Name of file to be written                */
  int           width,          /* Width of the image                        */
  int           height,         /* Height of the image                       */
  int           bx,             /* Boundary in X                             */
  int           by,             /* Boundary in Y                             */
  float         ***image,       /* Input image                               */
  int           bytes,          /* Byte format to use [1,2]                  */
  const char    *comment,       /* Comment lines                             */
  va_list       args            /* variable length argument list             */
)
{
  FILE         *image_file;
  va_list      args_copy;
  int          written;
  char         *filename_replaced;
  
  va_copy(args_copy,args);
  
  if (image == 0)
  {
    console_error("Could not write image %s - NULL pointer! Maybe there\n",
                  filename);
    console_text ("was a problem loading the image already?\n");
    return;
  }

  if ((bytes < 1) || (bytes > 2))
  {
    console_error("Unexpected format for file %s!\n", filename);
    return;
  }

  filename_replaced = (char*)calloc(1024,sizeof(char));
  written = vsnprintf (filename_replaced,1024,filename, args);
  if (!(image_file = fopen(filename_replaced, "w")))
  {
    console_error("Could not open file %s!\n", filename);
    return;
  }

  /* Write magic bytes                                                       */
  fprintf(image_file, "P6\n");

  /* Write (optional) comment lines                                          */
  if (comment)
  {
    int pos     = 0;
    int newline = 1;
    char *fullstring;
    char *comment_replaced;
	char *format;

	format = (char*)calloc(strlen(filename)+strlen(comment)+1,sizeof(char));
	strcpy(format,filename);
	strcat(format,comment);
	fullstring = (char*)calloc(2048, sizeof(char));
    vsnprintf (fullstring, 2048, format, args_copy);
	comment_replaced=fullstring+written;
	
    while (comment[pos] != '\0')
    {
      if (newline && (comment[pos] != '#'))
      {
        fprintf(image_file, "#%c", comment[pos]);
        newline = 0;
      }
      else
        fprintf(image_file, "%c", comment[pos]);

      if (comment[pos] == '\n')
        newline = 1;

      ++pos;
    }

    if (!newline)
      fprintf(image_file, "\n");
	
	free(fullstring);free(format);
  }

  /* Write dimensions                                                        */
  fprintf(image_file, "%d %d\n", width, height);

  /* Write depth                                                             */
  if (bytes == 1)
    fprintf(image_file, "255\n");
  else if (bytes == 2)
    fprintf(image_file, "65535\n");

  /* Read in the raw input                                                   */
  if (bytes == 1)
  {
    int i, j, c;

    for (j = by; j < height + by; ++j)
      for (i = bx; i < width + bx; ++i)
        for (c = 0; c < 3; ++c)
          fprintf(image_file, "%c", float_to_uchar(image[c][i][j]));
  }
  else if (bytes == 2)
  {
    int i, j, c;

    for (j = by; j < height + by; ++j)
      for (i = bx; i < width + bx; ++i)
        for (c = 0; c < 3; ++c)
          fprintf(image_file, "%2c",float_to_uchar(image[c][i][j]));
  }

  fclose(image_file);
  va_end (args_copy);
  printf("written %s\n",filename_replaced);
  free(filename_replaced);
}

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Write a one-byte or two-byte PPM file.                                    */
/*                                                                           */
/* Comment takes several comment lines, separated by '\n'. There is no need  */
/* to include markers such as '#', since they are automatically inserted     */
/* (where necessary)                                                         */
/* O.Demetz, June 2011: Extended functionality of printf-like features for   */
/* the filename and comment. Just give filename and comment with %s and %f   */
/* and such patterns, and append the corresponding variables like for printf.*/
/*****************************************************************************/
void write_ppm_image_comment
(
  const char    *filename,      /* Name of file to be written                */
  int           width,          /* Width of the image                        */
  int           height,         /* Height of the image                       */
  int           bx,             /* Boundary in X                             */
  int           by,             /* Boundary in Y                             */
  float         ***image,       /* Input image                               */
  int           bytes,          /* Byte format to use [1,2]                  */
  const char    *comment,       /* Comment lines                             */
  ...                           /* variable length argument list             */
)
{
	va_list args;
	va_start(args,comment);
	write_ppm_image_comment_started(filename,width,height,bx,by,image,bytes,
									comment,args);
	va_end(args);
}

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Write a one-byte or two-byte PPM file.                                    */
/*****************************************************************************/
void write_ppm_image
(
  const char    *filename,      /* Name of file to be written                */
  int           width,          /* Width of the image                        */
  int           height,         /* Height of the image                       */
  int           bx,             /* Boundary in X                             */
  int           by,             /* Boundary in Y                             */
  float         ***image,       /* Input image                               */
  int           bytes,          /* Byte format to use [1,2]                  */
  ...
)
{
  va_list args;
  va_start(args,bytes);
  write_ppm_image_comment_started(filename, width, height, bx, by, image, 
								  bytes, NULL, args);
  va_end(args);
}

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Free an image array allocated by read_pgm_image.                          */
/*****************************************************************************/
void free_pgm_image
(
  float         **image,        /* Input image                               */
  int           width,          /* Width of the image                        */
  int           height,         /* Height of the image                       */
  int           bx,             /* Boundary in X                             */
  int           by              /* Boundary in Y                             */
)
{
  if (image == 0)
  {
    console_error("Could not free image - NULL pointer! Maybe there\n");
    console_text ("was a problem loading the image already?\n");
    return;
  }

  free_multi(1,      2,       sizeof(float),
             width,  height,
             bx,     by,
             0,      0,
             &image);
}

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Free an image array allocated by read_ppm_image.                          */
/*****************************************************************************/
void free_ppm_image
(
  float         ***image,       /* Input image                               */
  int           width,          /* Width of the image                        */
  int           height,         /* Height of the image                       */
  int           bx,             /* Boundary in X                             */
  int           by              /* Boundary in Y                             */
)
{
  if (image == 0)
  {
    console_error("Could not free image - NULL pointer! Maybe there\n");
    console_text ("was a problem loading the image already?\n");
    return;
  }

  free_multi(1, 3,     sizeof(float),
             3, width, height,
             0, bx,    by,
             0, 0,     0,
             &image);
}

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Detect the type of a given image file.                                    */
/*****************************************************************************/
image_t detect_image_type
(
  const char    *filename       /* Name of file to be explored               */
)
{
  FILE         *image_file;
  char         *line = NULL;
  size_t       length = 0;

  if (!(image_file = fopen(filename, "r")))
  {
    console_error("Could not open file %s!\n", filename);
    return IMG_TYPE_INVALID;
  }

  /* Check for the magic bytes                                               */
  if (getline(&line, &length, image_file) >= 2)
  {
    if ((line[0] == 'P') && (line[1] == '5'))
    {
      fclose(image_file);
      return IMG_TYPE_PGM;
    }
    else if ((line[0] == 'P') && (line[1] == '6'))
    {
      fclose(image_file);
      return IMG_TYPE_PPM;
    }
  }

  fclose(image_file);
  return IMG_TYPE_INVALID;
}

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Duplicates a 2D image array.                                              */
/* RETURNS a pointer to the float-based image contents with bx,by wide       */
/*         boundaries (same initialisation as original image), or NULL on    */
/*         error.                                                            */
/*                                                                           */
/* The allocated memory needs to be freed by the user program afterwards,    */
/* either manually or by calling free_pgm_image.                             */
/*****************************************************************************/
float **duplicate_pgm_image
(
  float **in,   /* Input image                                               */
  int   nx,     /* Width                                                     */
  int   ny,     /* Height                                                    */
  int   bx,     /* Boundary in X                                             */
  int   by      /* Boundary in Y                                             */
)
{
  float **out;

  /* Allocate the image array                                                */
  if (!malloc_multi(1,  2,  sizeof(float),
                    nx, ny,
                    bx, by,
                    0,  0,
                    &out))
  {
    console_error("Cannot allocate memory to duplicate image!");
    return NULL;
  }

  memcpy(out[0], in[0], (nx + 2*bx) * (ny + 2*by) * sizeof(float));

  return out;
}

void read_barron_data
(
    char *filename, /* in : file name */
    float **u,     /* in : x-component of vector data */
    float **v,     /* in : y-component of vector data */
    int  nx,      /* in : size in x-direction */
    int  ny,      /* in : size in y-direction */
    int  bx,      /* in : boundary in x-direction */
    int  by       /* in : boundary in y-direction */      
)

/* reads barron file */

{
    FILE *file;    /* file pointer */
    float help;    /* tmp variable */
    float **uref;  /* tmp array */
    float **vref;  /* tmp array */
    int i,j;     /* loop variabeles */
    int nxref_and_offsetx; /* size in x-direction with crop offset */
    int nyref_and_offsety; /* size in y-direction with crop offset */
    int nxref_no_offsetx;  /* size in x-direction without crop offset */
    int nyref_no_offsety;  /* size in y-direction without crop offset */
    int offsetx; /* crop offset in x-direction */
    int offsety; /* crop offset in y-direction */


//     printf("\n Trying to read barron %s ...",filename);
    
    /* try to open file */
    file = fopen(filename,"r");
    
    /* if file not found */
    if (file==NULL)
    {
// 	printf("... FAILED");
	printf("\n\n PROGRAM ABORTED !!! \n\n");
	exit(0);
    }

    /* read header */
    fread (&help, sizeof(float), 1, file);
    nxref_and_offsetx  = (int) help;
    fread (&help, sizeof(float), 1, file);
    nyref_and_offsety  = (int) help;
    fread (&help, sizeof(float), 1, file);
    nxref_no_offsetx  = (int) help;
    fread (&help, sizeof(float), 1, file);
    nyref_no_offsety  = (int) help;
    fread (&help, sizeof(float), 1, file);
    offsetx = (int) help;
    fread (&help, sizeof(float), 1, file);
    offsety = (int) help;

    /* compare dimensions for consistency */
    if ((nx!=nxref_no_offsetx)||(ny!=nyref_no_offsety))
    {
// 	printf("... WRONG DIMENSIONS");
	printf("\n\n PROGRAM ABORTED !!! \n\n");
	exit(0);
    }

    /* allocate memory for tmp array */
    malloc_multi(2,2,sizeof(float),nxref_and_offsetx,nyref_and_offsety,0,0,0,0,&uref,&vref);

    /* read barron data */
    for (j = 0; j < nyref_and_offsety; j++)
	for (i = 0; i < nxref_and_offsetx; i++)
	{
	    fread(&help, sizeof(float), 1, file);	   
	    uref[i][j] = (float) help;	    
	    fread(&help, sizeof(float), 1, file);	   
	    vref[i][j] = (float) help;	   
	}    


    /* copy data without cropped border */
    for (i=bx; i<nx+bx; i++)
	for (j=by; j<ny+by; j++)
	{
	    u[i][j] = (float) uref[i-bx+offsetx][j-by+offsety];    
	    v[i][j] = (float) vref[i-bx+offsetx][j-by+offsety]; 
	}                        
    
    /* free memory for tmp array */
    free_multi(2,2,sizeof(float),nxref_and_offsetx,nyref_and_offsety,0,0,0,0,&uref,&vref);


    /* close file */
    fclose(file);

//     printf("... SUCCESS");
}

/* -------------------------------------------------------------------------- */
void write_barron_data
(
    char *filename, /* in : file name */
    float **u,     /* in : x-component of vector data */
    float **v,     /* in : y-component of vector data */
    int  nx,      /* in : size in x-direction */
    int  ny,      /* in : size in y-direction */
    int  bx,      /* in : boundary in x-direction */
    int  by       /* in : boundary in y-direction */       
)

/* writes barron file */

{
    FILE *file;   /* file pointer */
    float help;   /* tmp variable */
    int i,j;    /* loop variables */ 
    int offset; /* border size to crop (set fixed to 0) */ 

//     printf("\n Trying to write barron file %s ...",filename);
    
    /* try to open file */
    file = fopen(filename,"w");
    
    /* if file not found */
    if (file==NULL)
    {
// 	printf("... FAILED");
	printf("\n\n PROGRAM ABORTED !!! \n\n");
	exit(0);
    }
    
    /* write header */
    help = nx;
    fwrite (&help, 4, 1, file);
    help = ny;
    fwrite (&help, 4, 1, file);
    offset=0;
    help = nx - 2 * offset;
    fwrite (&help, 4, 1, file);
    help = ny - 2 * offset;
    fwrite (&help, 4, 1, file);
    help = offset;
    fwrite (&help, 4, 1, file);
    fwrite (&help, 4, 1, file);

    /* write data */
    for (j=by; j<ny+by; j++)
	for (i=bx; i<nx+bx; i++)
	{
	    help = (float)u[i][j];
	    fwrite (&help, 4, 1, file);
	    
	    help = (float)v[i][j];
	    fwrite (&help, 4, 1, file);
	}
       
    /* close file */
    fclose(file);
    
//     printf("... SUCCESS");
}

/* ------------------------------------------------------------------------- */

#endif

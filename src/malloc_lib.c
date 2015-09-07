#ifndef MALLOC_LIB
#define MALLOC_LIB
/*****************************************************************************/
/* --- malloc_lib ---------------------------------------------------------- */
/* Allocate complex structures for generic problems.                         */
/*                                                                           */
/* (C) 2008 - 2010 Pascal Gwosdek                                            */
/*                 Mathematical Image Analysis Group,                        */
/*                 Saarland University, Germany                              */
/*                                                                           */
/* Version 1.01 (2010-09-15)                                                 */
/*****************************************************************************/

#include <stdlib.h>
#include <stdarg.h>

#include "console_lib.c"


/* Prototypes for private routines (to suppress warnings)                    */
int  malloc_multi_create_hierarchy
     (int, int, int, int*, int*, void*, int, int);
void free_multi_destroy_hierarchy (void**, int, int);

/*****************************************************************************/
/* Some preprocessor stuff for the malloc_multi routine...                   */
/*****************************************************************************/
#ifndef NDEBUG
#define PG_WRONG_MALLOC_MULTI_CHECK                                           \
  if ((p[0] == 0) || (r[0] != 0))                                             \
  {                                                                           \
    console_error("Wrong call to malloc_multi!");                             \
    console_info ("Make sure you call this routine with pointers to the");    \
    console_text ("objects to be alloc'd, and remember the first");           \
    console_text ("restriction flag has to be zero anyway. Did you remember");\
    console_text ("to specify restriction flags at all?");                    \
  }
#else
#define PG_WRONG_MALLOC_MULTI_CHECK
#endif

/*****************************************************************************/
/* Create a pointer hierarchy over the given pointer, and relink it to the   */
/* new top level                                                             */
/*                                                                           */
/* RETURNS 1 on success and 0 on failure                                     */
/*****************************************************************************/
int malloc_multi_create_hierarchy
(
  int   number,         /* number of similar structures to allocate          */
  int   dimensions,     /* depth of each structure                           */
  int   remaining,      /* remaining number of elements for this subtree     */
  int   *n,             /* list of widths of the single stages               */
  int   *b,             /* list of bounds of the single stages               */
  void  *p,             /* list of pointers to the different structures      */
  int   firstcall,      /* are we in the 'innermost' dimension?              */
  int   elemsize        /* size of one element                               */
)
{
  char  ***pointer = (char***)p;
  char  **levelmem;
  int   levelsize, stepsize;
  int   i;

  /* Allocate enough pointer space for the now finest level                  */
  stepsize  = n[dimensions-1] + 2*b[dimensions-1];
  remaining = remaining / stepsize;

  if (!firstcall)
    stepsize *= sizeof(void*);
  else
    stepsize *= elemsize;

  levelsize = number * remaining;
  levelmem  = (char**)malloc(levelsize * sizeof(char*));

  if (!levelmem)
  {
    console_error("Cannot allocate memory hierarchy in dimensions=%d!",
                  dimensions);
    return 0;
  }

  /* Distribute the pointers onto the finer level                            */
  levelmem[0] = **pointer;
  for (i = 1; i < levelsize; ++i)
    levelmem[i] = &(levelmem[i-1][stepsize]);

  /* Relink the topmost layer                                                */
  for (i = 0; i < number; ++i)
    *pointer[i] = (char*)&(levelmem[remaining * i]);

  /* Recurse with the next coarser level                                     */
  if (dimensions > 2)
  {
    if (!malloc_multi_create_hierarchy(number, dimensions-1, remaining,
                                       n, b, p, 0, elemsize))
    {
      free(levelmem);
      return 0;
    }
  }

  return 1;
}

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Create a complex data structure of arbitrary elementary type              */
/*                                                                           */
/* SYNTAX:                                                                   */
/*         number_of_structures,                                             */
/*         dimensions,                                                       */
/*         size of one element,                                              */
/*         list of widths,                                                   */
/*         list of bounds,                                                   */
/*         list of scale flags,                                              */
/*         list of pointer addresses to store items to                       */
/*                                                                           */
/* EXAMPLE 1: 2 cubixes, 3x4x5, with 1-bounds                                */
/*                                                                           */
/*            malloc_multi(2, 3, sizeof(float), // 2 cubixes                 */
/*                         3, 4, 5,             // nx, ny, nz                */
/*                         1, 1, 1,             // bx, by, bz                */
/*                         0, 0, 0,             // no pyramimds              */
/*                         &cubix1, &cubix2);   // each of type float***     */
/*                                                                           */
/* EXAMPLE 2: 2 matrix pyramids, 4x4, depth 3, fully scaling, with 1-bounds  */
/*                                                                           */
/*            malloc_multi(2, 3, sizeof(int),   // 2 matrix pyramids         */
/*                         3, 4, 4,             // max_rec_depth, nx, ny     */
/*                         0, 1, 1,             // -, bx, by                 */
/*                         0, 1, 1,             // x and y are scaled        */
/*                         &mpyr1, &mpyr2);     // each of type int***       */
/*                                                                           */
/* NOTE: Make sure you add one to the pyramid depth if porting code from     */
/*       Andres' libs to this function!!                                     */
/*                                                                           */
/* RETURNS 1 on success and 0 on failure                                     */
/*****************************************************************************/
#define PG_GENERATE_MALLOC_MULTI(PG_MALLOC_TYPE,PG_MALLOC_CONTENT)            \
int PG_MALLOC_TYPE##_multi                                                    \
(                                                                             \
  int   number,         /* Number of structures to allocate                */ \
  int   dimensions,     /* Number of dimensions per structure              */ \
  int   elementsize,    /* Size of one element, use sizeof() as argument   */ \
  ...                   /* List of n's, b's, r's and the pointers          */ \
)                                                                             \
{                                                                             \
  va_list       arguments;            /* argument list                     */ \
  int           n   [dimensions];     /* list of widths                    */ \
  int           **npyr = NULL;        /* scaled version thereof (for pyr's)*/ \
  int           b   [dimensions];     /* list of bounds                    */ \
  int           r   [dimensions];     /* list of scale flags               */ \
  void          **p [number];         /* list of pointers to the objects   */ \
  char          *totalmem;            /* finest dimension (containing data)*/ \
  char          **levelmem;           /* coarsest dimension (contain. ptrs)*/ \
  int           i, j;                 /* loop counters                     */ \
  int           is_pyramid;           /* flag determining the pyramid case */ \
  int           totalsize, layersize; /* sizes for ~mem                    */ \
  char          **helppointer;        /* helper                            */ \
  int           *offset = NULL;       /* widths of the pyramid levels      */ \
                                                                              \
  va_start(arguments, elementsize);                                           \
                                                                              \
  /* Parse n list                                                          */ \
  for (i = 0; i < dimensions; ++i)                                            \
    n[i] = va_arg(arguments,int);                                             \
                                                                              \
  /* Parse b list                                                          */ \
  for (i = 0; i < dimensions; ++i)                                            \
    b[i] = va_arg(arguments,int);                                             \
                                                                              \
  /* Parse r list                                                          */ \
  is_pyramid = 0;                                                             \
  for (i = 0; i < dimensions; ++i)                                            \
  {                                                                           \
    r[i] = va_arg(arguments,int);                                             \
                                                                              \
    if (r[i])                                                                 \
      is_pyramid = 1;                                                         \
  }                                                                           \
                                                                              \
  /* Parse pointer list                                                    */ \
  for (i = 0; i < number; ++i)                                                \
    p[i] = va_arg(arguments,void**);                                          \
                                                                              \
  va_end(arguments);                                                          \
                                                                              \
  /* Check for typical wrong calls                                         */ \
  PG_WRONG_MALLOC_MULTI_CHECK                                                 \
                                                                              \
  /* Compute total amount of storage needed                                */ \
  if (is_pyramid)                                                             \
  {                                                                           \
    offset = (int*)malloc(n[0] * sizeof(int));                                \
                                                                              \
    if (!offset)                                                              \
    {                                                                         \
      console_error("Cannot allocate offset for pyramid!");                   \
      return 0;                                                               \
    }                                                                         \
                                                                              \
    npyr = (int**)malloc((n[0] + 1) * sizeof(int*));                          \
                                                                              \
    if (!npyr)                                                                \
    {                                                                         \
      console_error("Cannot allocate npyr for pyramid!");                     \
      free(offset);                                                           \
      return 0;                                                               \
    }                                                                         \
                                                                              \
    for (i = 0; i <= n[0]; ++i)                                               \
    {                                                                         \
      npyr[i] = (int*)malloc(dimensions * sizeof(int));                       \
                                                                              \
      if (!npyr[i])                                                           \
      {                                                                       \
        console_error("Cannot allocate npyr for pyramid!");                   \
                                                                              \
        for (j = 0; j < i; ++j)                                               \
          free(npyr[j]);                                                      \
                                                                              \
        free(npyr);                                                           \
        free(offset);                                                         \
                                                                              \
        return 0;                                                             \
      }                                                                       \
    }                                                                         \
                                                                              \
    totalsize = 0;                                                            \
    for (i = 1; i < dimensions; ++i)                                          \
      npyr[0][i] = n[i];                                                      \
                                                                              \
    for (i = 0; i < n[0]; ++i)                                                \
    {                                                                         \
      layersize = 1;                                                          \
      for (j = 1; j < dimensions; ++j)                                        \
      {                                                                       \
        layersize *= npyr[i][j] + 2 * b[j];                                   \
        if (r[j])                                                             \
          npyr[i+1][j] = npyr[i][j] / 2 + npyr[i][j] % 2;                     \
        else                                                                  \
          npyr[i+1][j] = n[j];                                                \
      }                                                                       \
      offset[i]  = layersize * elementsize;                                   \
      totalsize += layersize;                                                 \
    }                                                                         \
  }                                                                           \
  else                                                                        \
  {                                                                           \
    totalsize = 1;                                                            \
                                                                              \
    for (i = 0; i < dimensions; ++i)                                          \
      totalsize *= n[i] + 2 * b[i];                                           \
  }                                                                           \
                                                                              \
  /* Allocate that much storage                                            */ \
  totalmem = (char*)PG_MALLOC_TYPE(PG_MALLOC_CONTENT);                        \
                                                                              \
  if (!totalmem)                                                              \
  {                                                                           \
    console_error("Cannot allocate main memory block!");                      \
                                                                              \
    if (is_pyramid)                                                           \
    {                                                                         \
      for (j = 0; j <= n[0]; ++j)                                             \
        free(npyr[j]);                                                        \
                                                                              \
      free(npyr);                                                             \
      free(offset);                                                           \
    }                                                                         \
                                                                              \
    return 0;                                                                 \
  }                                                                           \
                                                                              \
  /* Now comes the tricky part. First, we link the first layer of each     */ \
  /* structure in a 'flat' manner to the memory, and will then lift it     */ \
  /* outside stepwise                                                      */ \
                                                                              \
  /* Create the pointer structure upon                                     */ \
  if (is_pyramid)                                                             \
  {                                                                           \
    /* This is the slightly harder case: Pyramids                          */ \
    levelmem = (char**)malloc(n[0] * number * sizeof(char**));                \
                                                                              \
    if (!levelmem)                                                            \
    {                                                                         \
      console_error("Cannot allocate level memory block!");                   \
                                                                              \
      free(totalmem);                                                         \
                                                                              \
      for (j = 0; j <= n[0]; ++j)                                             \
        free(npyr[j]);                                                        \
                                                                              \
      free(npyr);                                                             \
      free(offset);                                                           \
                                                                              \
      return 0;                                                               \
    }                                                                         \
                                                                              \
    /* The special alignment of the pyramids in memory (one after the      */ \
    /* other, i.e. all levels of one pyramid consecutively) forces us to   */ \
    /* really work with offsets instead of a fast iteration.               */ \
                                                                              \
    /* First link to the different layers                                  */ \
                                                                              \
    /* Topmost layer                                                       */ \
    for (i = 0; i < number; ++i)                                              \
    {                                                                         \
      levelmem[i * n[0]] = &(totalmem[i * totalsize * elementsize]);          \
      *p[i] = levelmem + i * n[0];                                            \
    }                                                                         \
                                                                              \
    /* All remaining layers                                                */ \
    for (i = 1; i < n[0]; ++i)                                                \
      for (j = 0; j < number; ++j)                                            \
        levelmem[j * n[0] + i] = &(levelmem[j * n[0] + i - 1]                 \
                                           [offset[i-1]]);                    \
                                                                              \
    /* Now that we have links to the various levels, allocate the pointer  */ \
    /* structures upon them                                                */ \
    for (i = 0; i < n[0]; ++i)                                                \
    {                                                                         \
      for (j = 0; j < number; ++j)                                            \
      {                                                                       \
        helppointer = &(levelmem[j * n[0] + i]);                              \
                                                                              \
        if (dimensions > 2)                                                   \
          if (!malloc_multi_create_hierarchy(1, dimensions - 1, offset[i],    \
                                             npyr[i] + 1, b + 1,              \
                                             &helppointer,                    \
                                             1, elementsize))                 \
          {                                                                   \
            console_error("Cannot allocate memory hierarchy!");               \
                                                                              \
            free(levelmem);                                                   \
            free(totalmem);                                                   \
                                                                              \
            for (j = 0; j <= n[0]; ++j)                                       \
              free(npyr[j]);                                                  \
                                                                              \
            free(npyr);                                                       \
            free(offset);                                                     \
                                                                              \
            return 0;                                                         \
          }                                                                   \
      }                                                                       \
    }                                                                         \
                                                                              \
    free(offset);                                                             \
    for (i = 0; i <= n[0]; ++i)                                               \
      free(npyr[i]);                                                          \
    free(npyr);                                                               \
  }                                                                           \
  else                                                                        \
  {                                                                           \
    /* This is the easy case: No pyramids, but 'cubical' objects           */ \
    for (i = 0; i < number; ++i)                                              \
      *p[i] = &(totalmem[i * totalsize * elementsize]);                       \
                                                                              \
    if (dimensions > 1)                                                       \
      if (!malloc_multi_create_hierarchy(number, dimensions, totalsize,       \
                                         n, b, p, 1, elementsize))            \
      {                                                                       \
        console_error("Cannot allocate memory hierarchy!");                   \
        free(totalmem);                                                       \
                                                                              \
        return 0;                                                             \
      }                                                                       \
  }                                                                           \
                                                                              \
  return 1;                                                                   \
}

/* Generate the malloc_multi and the calloc_multi function                   */
#define PG_COMMA ,
PG_GENERATE_MALLOC_MULTI(malloc, totalsize * number * elementsize)
PG_GENERATE_MALLOC_MULTI(calloc, elementsize PG_COMMA totalsize * number)

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Internally needed to destroy a complex data structure.                    */
/*****************************************************************************/
void free_multi_destroy_hierarchy
(
  void  **p,            /* list of pointers to the structures                */
  int   dimensions,     /* number of dimensions to dispose                   */
  int   destroy_finest  /* for pyramids, this must only be called once!      */
)
{
  if (dimensions > 0)
    free_multi_destroy_hierarchy((void**)*p, dimensions - 1, destroy_finest);

  if (destroy_finest || (dimensions > 0))
    free(p);
}

/* ------------------------------------------------------------------------- */

/*****************************************************************************/
/* Destroy complex data structure                                            */
/*                                                                           */
/* SYNTAX: *Exactly* as for malloc_multi.                                    */
/*****************************************************************************/
void free_multi
(
  int   number,         /* Number of structures to free                      */
  int   dimensions,     /* Number of dimensions per structure                */
  int   elementsize,    /* Size of one element, use sizeof() as argument     */
  ...                   /* List of the pointers                              */
)
{
  va_list       arguments;              /* argument list                     */
  int           n   [dimensions];       /* list of widths                    */
  int           b   [dimensions];       /* list of bounds                    */
  int           r   [dimensions];       /* list of scale flags               */
  void          ***p [number];          /* list of pointers to the objects   */
  int           i;                      /* loop counter                      */
  int           is_pyramid;             /* flag determining the pyramid case */
                                        
  va_start(arguments, elementsize);

  /* Parse n list                                                            */
  for (i = 0; i < dimensions; ++i)
    n[i] = va_arg(arguments,int);

  /* Parse b list                                                            */
  for (i = 0; i < dimensions; ++i)
    b[i] = va_arg(arguments,int);

  /* Parse r list                                                            */
  is_pyramid = 0;
  for (i = 0; i < dimensions; ++i)
  {
    r[i] = va_arg(arguments,int);

    if (r[i])
      is_pyramid = 1;
  }

  /* Parse pointer list                                                      */
  for (i = 0; i < number; ++i)
    p[i] = va_arg(arguments,void***);

  va_end(arguments);

  /* Check for typical wrong calls                                           */
  #ifndef NDEBUG
  if ((p[0] == 0) || (r[0] != 0))
  {
    console_error("Wrong call to free_multi!");
    console_info ("Make sure you call this routine with pointers to the");
    console_info ("objects to be freed, and remember the first restriction");
    console_info ("flag has to be zero anyway. Did you remember to specify");
    console_info ("restriction flags at all?");
  }
  #endif

  if (is_pyramid)
  {
    free_multi_destroy_hierarchy((void**)***p, dimensions - 2, 1);

    for (i = 1; i < n[0] * number; ++i)
      free_multi_destroy_hierarchy((void**)((**p)[i]), dimensions - 2, 0);

    free(**p);
  }   
  else
    free_multi_destroy_hierarchy(**p, dimensions - 1, 1);
}

/* ------------------------------------------------------------------------- */


#endif

/*!
 ************************************************************************
 * \file  memalloc.h
 *
 * \brief
 *    Memory allocation and free helper funtions
 *
 * \author
 *    Main contributors (see contributors.h for copyright, address and affiliation details)
 *     - Karsten Suehring
 *     - Alexis Michael Tourapis         <alexismt@ieee.org> 
 *     - Yuwen He                        <yhe@dolby.com>
 *
 ************************************************************************
 */

#ifndef _MEMALLOC_H_
#define _MEMALLOC_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "global.h"
#include "dpb.h"

extern int  get_mem2Dmp  (PicMotionParams ***array2D, int dim0, int dim1);

extern byte** new_mem2D(int dim0, int dim1);
extern int  get_mem2D(byte ***array2D, int dim0, int dim1);
extern int  get_mem3D(byte ****array3D, int dim0, int dim1, int dim2);
extern int  get_mem4D(byte *****array4D, int dim0, int dim1, int dim2, int dim3);

extern int** new_mem2Dint(int dim0, int dim1);
extern int  get_mem2Dint(int ***array2D, int dim0, int dim1);
extern int  get_mem3Dint(int ****array3D, int dim0, int dim1, int dim2);
extern int  get_mem4Dint(int *****array4D, int dim0, int dim1, int dim2, int dim3);

extern int  get_mem1Dpel(imgpel **array2D, int dim0);
extern int  get_mem2Dpel(imgpel ***array2D, int dim0, int dim1);
extern int  get_mem2Dpel_pad(imgpel ***array2D, int dim0, int dim1, int iPadY, int iPadX);
extern int  get_mem3Dpel    (imgpel ****array3D, int dim0, int dim1, int dim2);
extern int  get_mem3Dpel_pad(imgpel ****array3D, int dim0, int dim1, int dim2, int iPadY, int iPadX);

extern void free_mem2Dmp   (PicMotionParams    **array2D);

extern void free_mem2D     (byte      **array2D);
extern void free_mem3D     (byte     ***array3D);
extern void free_mem4D     (byte    ****array4D);

extern void free_mem2Dint  (int       **array2D);
extern void free_mem3Dint  (int      ***array3D);
extern void free_mem4Dint  (int     ****array4D);

extern void free_mem1Dpel    (imgpel     *array1D);
extern void free_mem2Dpel    (imgpel    **array2D);
extern void free_mem2Dpel_pad(imgpel **array2D, int iPadY, int iPadX);
extern void free_mem3Dpel    (imgpel   ***array3D);
extern void free_mem3Dpel_pad(imgpel ***array3D, int iDim12, int iPadY, int iPadX);

extern int  init_top_bot_planes(imgpel **imgFrame, int height, imgpel ***imgTopField, imgpel ***imgBotField);
extern void free_top_bot_planes(imgpel **imgTopField, imgpel **imgBotField);

extern void no_mem_exit(const char *where);


static inline void* mem_malloc(size_t nitems)
{
  void *d;
  if((d = malloc(nitems)) == NULL)
  {
    no_mem_exit("malloc failed.\n");
    return NULL;
  }
  return d;
}

/*!
 ************************************************************************
 * \brief
 *    allocate and set memory aligned at SSE_MEMORY_ALIGNMENT
 *
 ************************************************************************/
static inline void* mem_calloc(size_t nitems, size_t size)
{
  size_t padded_size = nitems * size; 
  void *d = mem_malloc(padded_size);
  memset(d, 0, (int)padded_size);
  return d;
}

static inline void mem_free(void *a)
{
  free_pointer(a);
}

#ifdef __cplusplus
}
#endif

#endif


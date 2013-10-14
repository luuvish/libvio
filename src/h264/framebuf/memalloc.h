#ifndef _MEMALLOC_H_
#define _MEMALLOC_H_


#include "global.h"
#include "dpb.h"
#include "picture.h"

extern int  get_mem2Dmp (pic_motion_params ***array2D, int dim0, int dim1);
extern void free_mem2Dmp(pic_motion_params **array2D);

extern int  get_mem1Dpel(imgpel **array2D, int dim0);
extern int  get_mem2Dpel(imgpel ***array2D, int dim0, int dim1);
extern int  get_mem2Dpel_pad(imgpel ***array2D, int dim0, int dim1, int iPadY, int iPadX);
extern int  get_mem3Dpel    (imgpel ****array3D, int dim0, int dim1, int dim2);

extern void free_mem1Dpel    (imgpel     *array1D);
extern void free_mem2Dpel    (imgpel    **array2D);
extern void free_mem2Dpel_pad(imgpel **array2D, int iPadY, int iPadX);
extern void free_mem3Dpel    (imgpel   ***array3D);

extern void no_mem_exit(const char *where);


#endif

#ifndef _MEMALLOC_H_
#define _MEMALLOC_H_


#include "global.h"
#include "dpb.h"
#include "frame_buffer.h"

extern int  get_mem2Dmp  (pic_motion_params ***array2D, int dim0, int dim1);

extern byte** new_mem2D(int dim0, int dim1);
extern int  get_mem2D(byte ***array2D, int dim0, int dim1);
extern int  get_mem3D(byte ****array3D, int dim0, int dim1, int dim2);

extern int** new_mem2Dint(int dim0, int dim1);
extern int  get_mem2Dint(int ***array2D, int dim0, int dim1);
extern int  get_mem3Dint(int ****array3D, int dim0, int dim1, int dim2);
extern int  get_mem4Dint(int *****array4D, int dim0, int dim1, int dim2, int dim3);

extern int  get_mem1Dpel(imgpel **array2D, int dim0);
extern int  get_mem2Dpel(imgpel ***array2D, int dim0, int dim1);
extern int  get_mem2Dpel_pad(imgpel ***array2D, int dim0, int dim1, int iPadY, int iPadX);
extern int  get_mem3Dpel    (imgpel ****array3D, int dim0, int dim1, int dim2);
extern int  get_mem3Dpel_pad(imgpel ****array3D, int dim0, int dim1, int dim2, int iPadY, int iPadX);

extern void free_mem2Dmp   (pic_motion_params    **array2D);

extern void free_mem2D     (byte      **array2D);
extern void free_mem3D     (byte     ***array3D);

extern void free_mem2Dint  (int       **array2D);
extern void free_mem3Dint  (int      ***array3D);
extern void free_mem4Dint  (int     ****array4D);

extern void free_mem1Dpel    (imgpel     *array1D);
extern void free_mem2Dpel    (imgpel    **array2D);
extern void free_mem2Dpel_pad(imgpel **array2D, int iPadY, int iPadX);
extern void free_mem3Dpel    (imgpel   ***array3D);
extern void free_mem3Dpel_pad(imgpel ***array3D, int iDim12, int iPadY, int iPadX);

extern void no_mem_exit(const char *where);


#endif

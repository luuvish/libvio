#include "global.h"
#include "memalloc.h"

/*!
 ************************************************************************
 * \brief
 *    Initialize 2-dimensional top and bottom field to point to the proper
 *    lines in frame
 *
 * \par Output:
 *    memory size in bytes
 ************************************************************************/
int init_top_bot_planes(imgpel **imgFrame, int dim0, imgpel ***imgTopField, imgpel ***imgBotField)
{
  int i;

  if((*imgTopField   = (imgpel**)mem_malloc((dim0>>1) * sizeof(imgpel*))) == NULL)
    no_mem_exit("init_top_bot_planes: imgTopField");

  if((*imgBotField   = (imgpel**)mem_malloc((dim0>>1) * sizeof(imgpel*))) == NULL)
    no_mem_exit("init_top_bot_planes: imgBotField");

  for(i = 0; i < (dim0>>1); i++)
  {
    (*imgTopField)[i] =  imgFrame[2 * i    ];
    (*imgBotField)[i] =  imgFrame[2 * i + 1];
  }

  return dim0 * sizeof(imgpel*);
}

 /*!
 ************************************************************************
 * \brief
 *    free 2-dimensional top and bottom fields without freeing target memory
 *
 * \par Output:
 *    memory size in bytes
 ************************************************************************/
void free_top_bot_planes(imgpel **imgTopField, imgpel **imgBotField)
{
  mem_free (imgTopField);
  mem_free (imgBotField);
}

/*!
 ************************************************************************
 * \brief
 *    Allocate 2D memory array -> PicMotionParams array2D[dim0][dim1]
 *
 * \par Output:
 *    memory size in bytes
 ************************************************************************/
int get_mem2Dmp(PicMotionParams ***array2D, int dim0, int dim1)
{
  int i;

  if((*array2D    = (PicMotionParams**)mem_malloc(dim0 *      sizeof(PicMotionParams*))) == NULL)
    no_mem_exit("get_mem2Dmp: array2D");
  if((*(*array2D) = (PicMotionParams* )mem_calloc(dim0 * dim1, sizeof(PicMotionParams ))) == NULL)
    no_mem_exit("get_mem2Dmp: array2D");

  for(i = 1 ; i < dim0; i++)
    (*array2D)[i] =  (*array2D)[i-1] + dim1;

  return dim0 * (sizeof(PicMotionParams*) + dim1 * sizeof(PicMotionParams));
}

/*!
 ************************************************************************
 * \brief
 *    free 2D memory array
 *    which was allocated with get_mem2Dmp()
 ************************************************************************
 */
void free_mem2Dmp(PicMotionParams **array2D)
{
  if (array2D)
  {
    if (*array2D)
      mem_free (*array2D);
    else 
      error ("free_mem2Dmp: trying to free unused memory",100);

    mem_free (array2D);
  } 
  else
  {
    error ("free_mem2Dmp: trying to free unused memory",100);
  }
}

/*!
 ************************************************************************
 * \brief
 *    Allocate 1D memory array -> imgpel array1D[dim0
 *
 * \par Output:
 *    memory size in bytes
 ************************************************************************/
int get_mem1Dpel(imgpel **array1D, int dim0)
{
  if((*array1D    = (imgpel*)mem_calloc(dim0,       sizeof(imgpel))) == NULL)
    no_mem_exit("get_mem1Dpel: array1D");

  return (sizeof(imgpel*) + dim0 * sizeof(imgpel));
}

/*!
 ************************************************************************
 * \brief
 *    Allocate 2D memory array -> imgpel array2D[dim0][dim1]
 *
 * \par Output:
 *    memory size in bytes
 ************************************************************************/
int get_mem2Dpel(imgpel ***array2D, int dim0, int dim1)
{
  int i;

  if((*array2D    = (imgpel**)mem_malloc(dim0 *        sizeof(imgpel*))) == NULL)
    no_mem_exit("get_mem2Dpel: array2D");
  if((*(*array2D) = (imgpel* )mem_malloc(dim0 * dim1 * sizeof(imgpel ))) == NULL)
    no_mem_exit("get_mem2Dpel: array2D");

  for(i = 1 ; i < dim0; i++)
  {
    (*array2D)[i] = (*array2D)[i-1] + dim1;
  }

  return dim0 * (sizeof(imgpel*) + dim1 * sizeof(imgpel));
}

int get_mem2Dpel_pad(imgpel ***array2D, int dim0, int dim1, int iPadY, int iPadX)
{
  int i;
  imgpel *curr = NULL;
  int iHeight, iWidth;
  
  iHeight = dim0+2*iPadY;
  iWidth = dim1+2*iPadX;
  if((*array2D    = (imgpel**)mem_malloc(iHeight*sizeof(imgpel*))) == NULL)
    no_mem_exit("get_mem2Dpel_pad: array2D");
  if((*(*array2D) = (imgpel* )mem_calloc(iHeight * iWidth, sizeof(imgpel ))) == NULL)
    no_mem_exit("get_mem2Dpel_pad: array2D");

  (*array2D)[0] += iPadX;
  curr = (*array2D)[0];
  for(i = 1 ; i < iHeight; i++)
  {
    curr += iWidth;
    (*array2D)[i] = curr;
  }
  (*array2D) = &((*array2D)[iPadY]);

  return iHeight * (sizeof(imgpel*) + iWidth * sizeof(imgpel));
}


/*!
 ************************************************************************
 * \brief
 *    Allocate 3D memory array -> imgpel array3D[dim0][dim1][dim2]
 *
 * \par Output:
 *    memory size in bytes
 ************************************************************************
 */
int get_mem3Dpel(imgpel ****array3D, int dim0, int dim1, int dim2)
{
  int i, mem_size = dim0 * sizeof(imgpel**);

  if(((*array3D) = (imgpel***)malloc(dim0 * sizeof(imgpel**))) == NULL)
    no_mem_exit("get_mem3Dpel: array3D");

  mem_size += get_mem2Dpel(*array3D, dim0 * dim1, dim2);

  for(i = 1; i < dim0; i++)
    (*array3D)[i] = (*array3D)[i - 1] + dim1;
  
  return mem_size;
}

int get_mem3Dpel_pad(imgpel ****array3D, int dim0, int dim1, int dim2, int iPadY, int iPadX)
{
  int i, mem_size = dim0 * sizeof(imgpel**);

  if(((*array3D) = (imgpel***)mem_malloc(dim0*sizeof(imgpel**))) == NULL)
    no_mem_exit("get_mem3Dpel_pad: array3D");

  for(i = 0; i < dim0; i++)
    mem_size += get_mem2Dpel_pad((*array3D)+i, dim1, dim2, iPadY, iPadX);
  
  return mem_size;
}


/*!
 ************************************************************************
 * \brief
 *    Allocate 4D memory array -> imgpel array4D[dim0][dim1][dim2][dim3]
 *
 * \par Output:
 *    memory size in bytes
 ************************************************************************
 */
int get_mem4Dpel(imgpel *****array4D, int dim0, int dim1, int dim2, int dim3)
{  
  int  i, mem_size = dim0 * sizeof(imgpel***);

  if(((*array4D) = (imgpel****)mem_malloc(dim0 * sizeof(imgpel***))) == NULL)
    no_mem_exit("get_mem4Dpel: array4D");

  mem_size += get_mem3Dpel(*array4D, dim0 * dim1, dim2, dim3);

  for(i = 1; i < dim0; i++)
    (*array4D)[i] = (*array4D)[i - 1] + dim1;

  return mem_size;
}


int get_mem4Dpel_pad(imgpel *****array4D, int dim0, int dim1, int dim2, int dim3, int iPadY, int iPadX)
{  
  int  i, mem_size = dim0 * sizeof(imgpel***);

  if(((*array4D) = (imgpel****)mem_malloc(dim0 * sizeof(imgpel***))) == NULL)
    no_mem_exit("get_mem4Dpel_pad: array4D");

  mem_size += get_mem3Dpel_pad(*array4D, dim0 * dim1, dim2, dim3, iPadY, iPadX);

  for(i = 1; i < dim0; i++)
    (*array4D)[i] = (*array4D)[i - 1] + dim1;

  return mem_size;
}

/*!
 ************************************************************************
 * \brief
 *    free 1D memory array
 *    which was allocated with get_mem1Dpel()
 ************************************************************************
 */
void free_mem1Dpel(imgpel *array1D)
{
  if (array1D)
  {
    mem_free (array1D);
  } 
  else
  {
    error ("free_mem1Dpel: trying to free unused memory",100);
  }
}

/*!
 ************************************************************************
 * \brief
 *    free 2D memory array
 *    which was allocated with get_mem2Dpel()
 ************************************************************************
 */
void free_mem2Dpel(imgpel **array2D)
{
  if (array2D)
  {
    if (*array2D)
      mem_free (*array2D);
    else 
     error ("free_mem2Dpel: trying to free unused memory",100);

    mem_free (array2D);
  } 
  else
  {
    error ("free_mem2Dpel: trying to free unused memory",100);
  }
}

void free_mem2Dpel_pad(imgpel **array2D, int iPadY, int iPadX)
{
  if (array2D)
  {
    if (*array2D)
    {
      mem_free (array2D[-iPadY]-iPadX);
    }
    else 
      error ("free_mem2Dpel_pad: trying to free unused memory",100);

    mem_free (&array2D[-iPadY]);
  } 
  else
  {
    error ("free_mem2Dpel_pad: trying to free unused memory",100);
  }
}

/*!
 ************************************************************************
 * \brief
 *    free 3D memory array
 *    which was allocated with get_mem3Dpel()
 ************************************************************************
 */
void free_mem3Dpel(imgpel ***array3D)
{
  if (array3D)
  {
    free_mem2Dpel(*array3D);
    mem_free (array3D);
  }
  else
  {
    error ("free_mem3Dpel: trying to free unused memory",100);
  }
}

void free_mem3Dpel_pad(imgpel ***array3D, int iDim12, int iPadY, int iPadX)
{
  if (array3D)
  {
    int i;
    for(i=0; i<iDim12; i++)
      if(array3D[i])
      {
        free_mem2Dpel_pad(array3D[i], iPadY, iPadX);
        array3D[i] = NULL;
      }
    mem_free (array3D);
  }
  else
  {
    error ("free_mem3Dpel_pad: trying to free unused memory",100);
  }
  
}

/*!
 ************************************************************************
 * \brief
 *    Create 2D memory array -> byte array2D[dim0][dim1]
 *
 * \par Output:
 *    byte type array of size dim0 * dim1
 ************************************************************************
 */
byte** new_mem2D(int dim0, int dim1)
{
  int i;
  byte **array2D;

  if((  array2D  = (byte**)mem_malloc(dim0 *      sizeof(byte*))) == NULL)
    no_mem_exit("get_mem2D: array2D");
  if((*(array2D) = (byte* )mem_calloc(dim0 * dim1,sizeof(byte ))) == NULL)
    no_mem_exit("get_mem2D: array2D");

  for(i = 1; i < dim0; i++)
    array2D[i] = array2D[i-1] + dim1;

  return (array2D);
}

/*!
 ************************************************************************
 * \brief
 *    Allocate 2D memory array -> unsigned char array2D[dim0][dim1]
 *
 * \par Output:
 *    memory size in bytes
 ************************************************************************/
int get_mem2D(byte ***array2D, int dim0, int dim1)
{
  int i;

  if((  *array2D  = (byte**)mem_malloc(dim0 *      sizeof(byte*))) == NULL)
    no_mem_exit("get_mem2D: array2D");
  if((*(*array2D) = (byte* )mem_calloc(dim0 * dim1,sizeof(byte ))) == NULL)
    no_mem_exit("get_mem2D: array2D");

  for(i = 1; i < dim0; i++)
    (*array2D)[i] = (*array2D)[i-1] + dim1;

  return dim0 * (sizeof(byte*) + dim1 * sizeof(byte));
}


/*!
 ************************************************************************
 * \brief
 *    Create 2D memory array -> int array2D[dim0][dim1]
 *
 * \par Output:
 *    int type array of size dim0 * dim1
 ************************************************************************
 */
int** new_mem2Dint(int dim0, int dim1)
{
  int i;
  int **array2D;

  if((array2D    = (int**)mem_malloc(dim0 *       sizeof(int*))) == NULL)
    no_mem_exit("get_mem2Dint: array2D");
  if((*(array2D) = (int* )mem_calloc(dim0 * dim1, sizeof(int ))) == NULL)
    no_mem_exit("get_mem2Dint: array2D");

  for(i = 1 ; i < dim0; i++)
    (array2D)[i] =  (array2D)[i-1] + dim1;

  return (array2D);
}

/*!
 ************************************************************************
 * \brief
 *    Allocate 2D memory array -> int array2D[dim0][dim1]
 *
 * \par Output:
 *    memory size in bytes
 ************************************************************************
 */
int get_mem2Dint(int ***array2D, int dim0, int dim1)
{
  int i;

  if((*array2D    = (int**)mem_malloc(dim0 *       sizeof(int*))) == NULL)
    no_mem_exit("get_mem2Dint: array2D");
  if((*(*array2D) = (int* )mem_calloc(dim0 * dim1, sizeof(int ))) == NULL)
    no_mem_exit("get_mem2Dint: array2D");

  for(i = 1 ; i < dim0; i++)
    (*array2D)[i] =  (*array2D)[i-1] + dim1;

  return dim0 * (sizeof(int*) + dim1 * sizeof(int));
}

/*!
 ************************************************************************
 * \brief
 *    Allocate 3D memory array -> unsigned char array3D[dim0][dim1][dim2]
 *
 * \par Output:
 *    memory size in bytes
 ************************************************************************
 */
int get_mem3D(byte ****array3D, int dim0, int dim1, int dim2)
{
  int  i, mem_size = dim0 * sizeof(byte**);

  if(((*array3D) = (byte***)mem_malloc(dim0 * sizeof(byte**))) == NULL)
    no_mem_exit("get_mem3D: array3D");

  mem_size += get_mem2D(*array3D, dim0 * dim1, dim2);

  for(i = 1; i < dim0; i++)
    (*array3D)[i] =  (*array3D)[i-1] + dim1;

  return mem_size;
}

/*!
 ************************************************************************
 * \brief
 *    Allocate 4D memory array -> unsigned char array4D[dim0][dim1][dim2][dim3]
 *
 * \par Output:
 *    memory size in bytes
 ************************************************************************
 */
int get_mem4D(byte *****array4D, int dim0, int dim1, int dim2, int dim3)
{
  int  i, mem_size = dim0 * sizeof(byte***);

  if(((*array4D) = (byte****)mem_malloc(dim0 * sizeof(byte***))) == NULL)
    no_mem_exit("get_mem4D: array4D");

  mem_size += get_mem3D(*array4D, dim0 * dim1, dim2, dim3);

  for(i = 1; i < dim0; i++)
    (*array4D)[i] =  (*array4D)[i-1] + dim1;

  return mem_size;
}

/*!
 ************************************************************************
 * \brief
 *    Allocate 3D memory array -> int array3D[dim0][dim1][dim2]
 *
 * \par Output:
 *    memory size in bytes
 ************************************************************************
 */
int get_mem3Dint(int ****array3D, int dim0, int dim1, int dim2)
{
  int  i, mem_size = dim0 * sizeof(int**);

  if(((*array3D) = (int***)mem_malloc(dim0 * sizeof(int**))) == NULL)
    no_mem_exit("get_mem3Dint: array3D");

  mem_size += get_mem2Dint(*array3D, dim0 * dim1, dim2);

  for(i = 1; i < dim0; i++)
    (*array3D)[i] =  (*array3D)[i-1] + dim1;

  return mem_size;
}

/*!
 ************************************************************************
 * \brief
 *    Allocate 4D memory array -> int array4D[dim0][dim1][dim2][dim3]
 *
 * \par Output:
 *    memory size in bytes
 ************************************************************************
 */
int get_mem4Dint(int *****array4D, int dim0, int dim1, int dim2, int dim3)
{
  int  i, mem_size = dim0 * sizeof(int***);

  if(((*array4D) = (int****)mem_malloc(dim0 * sizeof(int***))) == NULL)
    no_mem_exit("get_mem4Dint: array4D");

  mem_size += get_mem3Dint(*array4D, dim0 * dim1, dim2, dim3);

  for(i = 1; i < dim0; i++)
    (*array4D)[i] =  (*array4D)[i-1] + dim1;

  return mem_size;
}


/*!
 ************************************************************************
 * \brief
 *    free 2D memory array
 *    which was allocated with get_mem2D()
 ************************************************************************
 */
void free_mem2D(byte **array2D)
{
  if (array2D)
  {
    if (*array2D)
      mem_free (*array2D);
    else 
      error ("free_mem2D: trying to free unused memory",100);

    mem_free (array2D);
  } 
  else
  {
    error ("free_mem2D: trying to free unused memory",100);
  }
}

/*!
 ************************************************************************
 * \brief
 *    free 2D memory array
 *    which was allocated with get_mem2Dint()
 ************************************************************************
 */
void free_mem2Dint(int **array2D)
{
  if (array2D)
  {
    if (*array2D)
      mem_free (*array2D);
    else 
      error ("free_mem2Dint: trying to free unused memory",100);

    mem_free (array2D);
  } 
  else
  {
    error ("free_mem2Dint: trying to free unused memory",100);
  }
}

/*!
 ************************************************************************
 * \brief
 *    free 3D memory array
 *    which was allocated with get_mem3D()
 ************************************************************************
 */
void free_mem3D(byte ***array3D)
{
  if (array3D)
  {
   free_mem2D(*array3D);
   mem_free (array3D);
  } 
  else
  {
    error ("free_mem3D: trying to free unused memory",100);
  }
}

/*!
 ************************************************************************
 * \brief
 *    free 4D memory array
 *    which was allocated with get_mem3D()
 ************************************************************************
 */
void free_mem4D(byte ****array4D)
{
  if (array4D)
  {
   free_mem3D(*array4D);
   mem_free (array4D);
  } 
  else
  {
    error ("free_mem4D: trying to free unused memory",100);
  }
}

/*!
 ************************************************************************
 * \brief
 *    free 3D memory array
 *    which was allocated with get_mem3Dint()
 ************************************************************************
 */
void free_mem3Dint(int ***array3D)
{
  if (array3D)
  {
   free_mem2Dint(*array3D);
   mem_free (array3D);
  } 
  else
  {
    error ("free_mem3Dint: trying to free unused memory",100);
  }
}

/*!
 ************************************************************************
 * \brief
 *    free 4D memory array
 *    which was allocated with get_mem4Dint()
 ************************************************************************
 */
void free_mem4Dint(int ****array4D)
{
  if (array4D)
  {
    free_mem3Dint( *array4D);
    mem_free (array4D);
  } else
  {
    error ("free_mem4Dint: trying to free unused memory",100);
  }
}

/*!
 ************************************************************************
 * \brief
 *    Exit program if memory allocation failed (using error())
 * \param where
 *    string indicating which memory allocation failed
 ************************************************************************
 */
void no_mem_exit(const char *where)
{
   snprintf(errortext, ET_SIZE, "Could not allocate memory: %s",where);
   error (errortext, 100);
}

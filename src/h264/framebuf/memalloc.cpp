#include "global.h"
#include "memalloc.h"


void no_mem_exit(const char *where)
{
   snprintf(errortext, ET_SIZE, "Could not allocate memory: %s",where);
   error (errortext, 100);
}


static inline void* mem_malloc(size_t nitems)
{
    void *d;
    if ((d = malloc(nitems)) == NULL) {
        no_mem_exit("malloc failed.\n");
        return NULL;
    }
    return d;
}

static inline void mem_free(void *pointer)
{
    if (pointer != NULL) {
        free(pointer);
        pointer = NULL;
    }
}

static inline void* mem_calloc(size_t nitems, size_t size)
{
    size_t padded_size = nitems * size; 
    void *d = mem_malloc(padded_size);
    memset(d, 0, (int)padded_size);
    return d;
}


int get_mem2Dmp(pic_motion_params ***array2D, int dim0, int dim1)
{
    if ((*array2D    = (pic_motion_params**)mem_malloc(dim0 *      sizeof(pic_motion_params*))) == NULL)
        no_mem_exit("get_mem2Dmp: array2D");
    if((*(*array2D) = (pic_motion_params* )mem_calloc(dim0 * dim1, sizeof(pic_motion_params ))) == NULL)
        no_mem_exit("get_mem2Dmp: array2D");

    for (int i = 1; i < dim0; ++i)
        (*array2D)[i] = (*array2D)[i - 1] + dim1;

    return dim0 * (sizeof(pic_motion_params*) + dim1 * sizeof(pic_motion_params));
}

void free_mem2Dmp(pic_motion_params** array2D)
{
    if (array2D) {
        if (*array2D)
            mem_free(*array2D);
        else
            error("free_mem2Dmp: trying to free unused memory", 100);
        mem_free(array2D);
    }
}

int get_mem1Dpel(px_t** array1D, int dim0)
{
    if ((*array1D    = (px_t*)mem_calloc(dim0,       sizeof(px_t))) == NULL)
        no_mem_exit("get_mem1Dpel: array1D");

    return (sizeof(px_t*) + dim0 * sizeof(px_t));
}

int get_mem2Dpel(px_t*** array2D, int dim0, int dim1)
{
    if ((*array2D    = (px_t**)mem_malloc(dim0 *        sizeof(px_t*))) == NULL)
        no_mem_exit("get_mem2Dpel: array2D");
    if ((*(*array2D) = (px_t* )mem_malloc(dim0 * dim1 * sizeof(px_t ))) == NULL)
        no_mem_exit("get_mem2Dpel: array2D");

    for (int i = 1; i < dim0; ++i)
        (*array2D)[i] = (*array2D)[i - 1] + dim1;

    return dim0 * (sizeof(px_t*) + dim1 * sizeof(px_t));
}

int get_mem2Dpel_pad(px_t*** array2D, int dim0, int dim1, int iPadY, int iPadX)
{
    int iHeight = dim0 + 2 * iPadY;
    int iWidth  = dim1 + 2 * iPadX;
    if ((*array2D    = (px_t**)mem_malloc(iHeight * sizeof(px_t*))) == NULL)
        no_mem_exit("get_mem2Dpel_pad: array2D");
    if ((*(*array2D) = (px_t* )mem_calloc(iHeight * iWidth, sizeof(px_t ))) == NULL)
        no_mem_exit("get_mem2Dpel_pad: array2D");

    (*array2D)[0] += iPadX;
    px_t* curr = (*array2D)[0];
    for (int i = 1; i < iHeight; ++i) {
        curr += iWidth;
        (*array2D)[i] = curr;
    }
    (*array2D) = &((*array2D)[iPadY]);

    return iHeight * (sizeof(px_t*) + iWidth * sizeof(px_t));
}

int get_mem3Dpel(px_t**** array3D, int dim0, int dim1, int dim2)
{
    int mem_size = dim0 * sizeof(px_t**);

    if (((*array3D) = (px_t***)malloc(dim0 * sizeof(px_t**))) == NULL)
        no_mem_exit("get_mem3Dpel: array3D");

    mem_size += get_mem2Dpel(*array3D, dim0 * dim1, dim2);

    for (int i = 1; i < dim0; ++i)
        (*array3D)[i] = (*array3D)[i - 1] + dim1;
  
    return mem_size;
}


void free_mem1Dpel(px_t* array1D)
{
    if (array1D)
        mem_free(array1D);
}

void free_mem2Dpel(px_t** array2D)
{
    if (array2D) {
        if (*array2D)
            mem_free(*array2D);
        else
            error("free_mem2Dpel: trying to free unused memory", 100);
        mem_free(array2D);
    }
}

void free_mem2Dpel_pad(px_t** array2D, int iPadY, int iPadX)
{
    if (array2D) {
        if (*array2D)
            mem_free(array2D[-iPadY] - iPadX);
        else
            error("free_mem2Dpel_pad: trying to free unused memory", 100);
        mem_free(&array2D[-iPadY]);
    }
}

void free_mem3Dpel(px_t*** array3D)
{
    if (array3D) {
        free_mem2Dpel(*array3D);
        mem_free(array3D);
    }
}

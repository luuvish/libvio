//
//  ocl_trquant.c
//  HM
//
//  Created by Injo Hwang on 12. 9. 29..
//
//

#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <OpenCL/opencl.h>

int              device_index = 0;
cl_device_type   device_type = CL_DEVICE_TYPE_GPU;
cl_device_id     device;
cl_context       context;
cl_command_queue queue;
char *           source;
cl_program       program;
cl_kernel        kernel;


void read_kernel(const char *name)
{
    int fd;
    unsigned length;
    struct stat status;
    int ret;

    fd = open(name, O_RDONLY);
    if (fd == -1) {
        printf("Error opening file %s\n", name);
        exit(1);
    }
    ret = fstat(fd, &status);
    if (ret) {
        printf("Error reading status for file %s\n", name);
        exit(1);
    }
    length = status.st_size;
    
    source = (char *)calloc(length+1, sizeof(char));
    ret = read(fd, source, length);
    if (!ret) {
        printf("Error reading from file %s\n", name);
        exit(1);
    }
    
    close(fd);
}

void init_opencl()
{
    cl_int err;
    cl_uint num_devices;

    // How many devices of the type requested are in the system?
    clGetDeviceIDs(NULL, device_type, 0, NULL, &num_devices);

    // Make sure the requested index is within bounds.  Otherwise, correct it.
    if (device_index < 0 || device_index > num_devices - 1) {
        printf("Requested index (%d) is out of range. Using 0.\n", device_index);
        device_index = 0;
    }

    // Grab the requested device.
    cl_device_id all_devices[num_devices];
    err = clGetDeviceIDs(NULL, CL_DEVICE_TYPE_GPU, num_devices, all_devices, NULL);
    device = all_devices[device_index];
    
    // Dump the device.
    char name[128];
    clGetDeviceInfo(device, CL_DEVICE_NAME, 128*sizeof(char), name, NULL);
    printf("Using OpenCL device: %s\n", name);
    
    // Create an OpenCL context using this compute device.
    context = clCreateContext(NULL, 1, &device, NULL, NULL, &err);
    if (!context) {
        printf("Error: Failed to create a compute context!\n");
        exit(1);
    }
    
    // Create a command queue on this device, since we want to use if for
    // running our CL program.
    queue = clCreateCommandQueue(context, device, 0, &err);
    if (!queue) {
        printf("Error: Failed to create a command queue!\n");
        exit(1);
    }

    read_kernel("/Users/injoh/Documents/workspace/video/hm-8.0.bak/source/Lib/TLibCommon/ocl_trquant.cl");

    program = clCreateProgramWithSource(context, 1, (const char **)&source, NULL, &err);
    if (!program) {
        printf("Error: Failed to create a compute program!\n");
        exit(1);
    }

    err = clBuildProgram(program, 1, &device, NULL, NULL, NULL);
    if (err != CL_SUCCESS) {
        size_t len;
        char buffer[2048];
        printf("Error: Failed to build program executable!\n");
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        printf("%s\n", buffer);
        exit(1);
    }
    
    kernel = clCreateKernel(program, "xITrMxN", &err);
    if (!kernel || err != CL_SUCCESS) {
        printf("Error: Failed to create a compute kernel!\n");
        exit(1);
    }
}

void close_opencl()
{
    clReleaseKernel(kernel);
    clReleaseProgram(program);
    clReleaseCommandQueue(queue);
    clReleaseContext(context);
}

void ocl_trquant(short *coeff, short *block,
                 const int iWidth, const int iHeight,
                 const unsigned int uiMode)
{
    size_t global = 1;

    cl_int err;
    cl_mem input = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                  sizeof(short) * (64*64), coeff, &err);
    cl_mem output = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
                                   sizeof(short) * (64*64), NULL, &err);

    if (!input || !output) {
        printf("Error: Failed to allocate device memory!\n");
        exit(1);
    }

    err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &input);
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &output);
    err |= clSetKernelArg(kernel, 2, sizeof(int), &iWidth);
    err |= clSetKernelArg(kernel, 3, sizeof(int), &iHeight);
    err |= clSetKernelArg(kernel, 4, sizeof(unsigned int), &uiMode);
    if (err != CL_SUCCESS) {
        printf("Error: Failed to set kernel arguments! %d\n", err);
        exit(1);
    }

    err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, &global, NULL, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        printf("Error: Failed to read output array! %d\n", err);
        exit(1);
    }

    err = clEnqueueReadBuffer(queue, output, CL_FALSE, 0, sizeof(short) * (64*64), block, 0, NULL, NULL);
    if (err != CL_SUCCESS) {
        printf("Error: Failed to read output array! %d\n", err);
        exit(1);
    }
    
    clFinish(queue);

    clReleaseMemObject(input);
    clReleaseMemObject(output);
}

#ifndef _WIN32_H_
#define _WIN32_H_


#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <time.h>
#include <stdint.h>
#if defined(OPENMP)
#include <omp.h>
#endif

#define TIME_T   struct timeval
typedef int64_t int64;

extern void  gettime(TIME_T* time);
extern void  init_time(void);
extern int64 timediff(TIME_T* start, TIME_T* end);
extern int64 timenorm(int64 cur_time);


#endif

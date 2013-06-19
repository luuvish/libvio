
/*!
 ************************************************************************
 *  \file
 *     win32.h
 *
 *  \brief
 *     win32 definitions for H.264 codec.
 *
 *  \author
 *
 ************************************************************************
 */
#ifndef _WIN32_H_
#define _WIN32_H_

#ifdef __cplusplus
extern "C" {
#endif

# include <fcntl.h>
# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <assert.h>
# include <unistd.h>
# include <sys/time.h>
# include <sys/stat.h>
# include <time.h>
# include <stdint.h>
#if defined(OPENMP)
# include <omp.h>
#endif

#define TIME_T   struct timeval
#define tell(fd) lseek(fd, 0, SEEK_CUR)
#define OPENFLAGS_WRITE O_WRONLY|O_CREAT|O_TRUNC
#define OPENFLAGS_READ  O_RDONLY
#define OPEN_PERMISSIONS S_IRUSR | S_IWUSR

typedef int64_t int64;
typedef uint64_t uint64;
# define FORMAT_OFF_T "lld"

extern void   gettime(TIME_T* time);
extern void   init_time(void);
extern int64 timediff(TIME_T* start, TIME_T* end);
extern int64 timenorm(int64 cur_time);

#ifdef __cplusplus
}
#endif

#endif

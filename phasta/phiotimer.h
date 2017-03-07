#ifndef PHIOTIMER_H
#define PHIOTIMER_H

/** \file phiotimer.h
    \brief timers for reading and writing phasta files
*/

#ifdef __cplusplus
extern "C" {
#endif

#define CHEFIO_READTIME(cmd,bytes) {\
    chefioTime t0,t1;\
    chefio_time(&t0);\
    cmd\
    chefio_time(&t1);\
    const size_t time = chefio_time_diff(&t0,&t1);\
    chefio_addReadTime(time);\
    chefio_addReadBytes(bytes);\
}

#define CHEFIO_WRITETIME(cmd,bytes) {\
    chefioTime t0,t1;\
    chefio_time(&t0);\
    cmd\
    chefio_time(&t1);\
    const size_t time = chefio_time_diff(&t0,&t1);\
    chefio_addWriteTime(time);\
    chefio_addWriteBytes(bytes);\
}

/* \brief constants to identify the different chef files */
enum chefio_file { CHEF_GEOMBC, CHEF_RESTART, NUM_CHEF_FILES};

#ifdef __INTEL_COMPILER
typedef size_t chefioTime;
#else
typedef struct timespec chefioTime;
#endif
/* \brief get the current time */
void chefio_time(chefioTime* t);
/* \brief compute the time difference end-start */
size_t chefio_time_diff(chefioTime* start, chefioTime* end);
/* \brief accumulate bytes read */
void chefio_addReadBytes(size_t b);
/* \brief accumulate bytes written */
void chefio_addWriteBytes(size_t b);
/* \brief accumulate time reading */
void chefio_addReadTime(size_t t);
/* \brief accumulate time writing */
void chefio_addWriteTime(size_t t);
/* \brief accumulate time opening */
void chefio_addOpenTime(size_t t);
/* \brief accumulate time closing */
void chefio_addCloseTime(size_t t);
/* \brief initialize the counters and timers */
void chefio_initStats();
/* \brief print io information */
void chefio_printStats();
/* \brief set the current file to record counters and timers for
 * \detail see chefio_file enum */
void chefio_setfile(int f);

#ifdef __cplusplus
}
#endif

#endif

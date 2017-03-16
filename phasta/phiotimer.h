#ifndef PHIOTIMER_H
#define PHIOTIMER_H

/** \file phiotimer.h
    \brief timers for reading and writing phasta files
*/

#ifdef __cplusplus
extern "C" {
#endif

#define PHASTAIO_READTIME(cmd,bytes) {\
    phastaioTime t0,t1;\
    phastaio_time(&t0);\
    cmd\
    phastaio_time(&t1);\
    const size_t time = phastaio_time_diff(&t0,&t1);\
    phastaio_addReadTime(time);\
    phastaio_addReadBytes(bytes);\
}

#define PHASTAIO_WRITETIME(cmd,bytes) {\
    phastaioTime t0,t1;\
    phastaio_time(&t0);\
    cmd\
    phastaio_time(&t1);\
    const size_t time = phastaio_time_diff(&t0,&t1);\
    phastaio_addWriteTime(time);\
    phastaio_addWriteBytes(bytes);\
}

#define PHASTAIO_OPENTIME(cmd) {\
    phastaioTime t0,t1;\
    phastaio_time(&t0);\
    cmd\
    phastaio_time(&t1);\
    const size_t time = phastaio_time_diff(&t0,&t1);\
    phastaio_addOpenTime(time);\
}

#define PHASTAIO_CLOSETIME(cmd) {\
    phastaioTime t0,t1;\
    phastaio_time(&t0);\
    cmd\
    phastaio_time(&t1);\
    const size_t time = phastaio_time_diff(&t0,&t1);\
    phastaio_addCloseTime(time);\
}

/* \brief constants to identify the different phasta and chef files */
enum phastaio_file {
  GEOMBC_READ,
  GEOMBC_WRITE,
  RESTART_READ,
  RESTART_WRITE,
  NUM_PHASTAIO_MODES
};

#ifdef __INTEL_COMPILER
typedef size_t phastaioTime;
#else
#include <time.h>
typedef struct timespec phastaioTime;
#endif
/* \brief get the current time */
void phastaio_time(phastaioTime* t);
/* \brief compute the time difference end-start */
size_t phastaio_time_diff(phastaioTime* start, phastaioTime* end);
/* \brief accumulate bytes read */
void phastaio_addReadBytes(size_t b);
/* \brief accumulate bytes written */
void phastaio_addWriteBytes(size_t b);
/* \brief accumulate time reading */
void phastaio_addReadTime(size_t t);
/* \brief accumulate time writing */
void phastaio_addWriteTime(size_t t);
/* \brief accumulate time opening */
void phastaio_addOpenTime(size_t t);
/* \brief accumulate time closing */
void phastaio_addCloseTime(size_t t);
/* \brief initialize the counters and timers */
void phastaio_initStats();
/* \brief print io information */
void phastaio_printStats();
/* \brief set the current file to record counters and timers for
 * \detail see phastaio_file enum */
void phastaio_setfile(int f);

#ifdef __cplusplus
}
#endif

#endif

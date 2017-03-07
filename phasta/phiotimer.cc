#include <stdio.h>
#include <cassert>
#include <phiotimer.h>
#include <PCU.h>
#include <inttypes.h> /* PRIu64 */
#include <time.h> /* clock_gettime */
#include <unistd.h> /* usleep */

#define BILLION 1000L*1000L*1000L
#define MILLION 1000L*1000L

struct chefio_stats {
  size_t cpus;
  size_t readTime[NUM_CHEF_FILES];
  size_t writeTime[NUM_CHEF_FILES];
  size_t readBytes[NUM_CHEF_FILES];
  size_t writeBytes[NUM_CHEF_FILES];
  size_t reads[NUM_CHEF_FILES];
  size_t writes[NUM_CHEF_FILES];
  size_t openTime[NUM_CHEF_FILES];
  size_t closeTime[NUM_CHEF_FILES];
  size_t opens[NUM_CHEF_FILES];
  size_t closes[NUM_CHEF_FILES];
  int fileIdx;
};
static struct chefio_stats chefio_global_stats;

#ifdef __INTEL_COMPILER
/* return the cycle count */
void chefio_time(chefioTime* t) {
  *t = _rdtsc(); //intel intrinsic
}
/* determine the reference clock frequency */
static size_t chefio_getCyclesPerMicroSec() {
  const size_t usec = 5*MILLION;
  size_t cpus, cycles;
  chefioTime t0, t1;
  chefio_time(&t0);
  /* Testing on Theta indicates that 5s is long enough
   * to get a stable value for the reference frequency.
   */
  usleep(usec);
  chefio_time(&t1);
  cycles = t1 - t0;
  cpus = ((double)cycles)/(usec);
  if(!PCU_Comm_Self())
    fprintf(stderr, "cycles %" PRIu64 " us %" PRIu64 " cycles per micro second %" PRIu64"\n", cycles, usec, cpus);
  return cpus;
}
/*return elapsed time in micro seconds*/
size_t chefio_time_diff(chefioTime* start, chefioTime* end) {
  size_t cycles = *end - *start;
  size_t us = ((double)cycles)/chefio_global_stats.cpus;
  return us;
}
#else
void chefio_time(chefioTime* t) {
  int err;
  err = clock_gettime(CLOCK_MONOTONIC,t);
  assert(!err);
}
/*return elapsed time in micro seconds*/
size_t chefio_time_diff(chefioTime* start, chefioTime* end) {
  assert(sizeof(size_t)==8);
  size_t elapsed = 0;
  chefioTime diff;
  if ((end->tv_nsec-start->tv_nsec)<0) {
    diff.tv_sec = end->tv_sec-start->tv_sec-1;
    diff.tv_nsec = BILLION+end->tv_nsec-start->tv_nsec;
  } else {
    diff.tv_sec = end->tv_sec-start->tv_sec;
    diff.tv_nsec = end->tv_nsec-start->tv_nsec;
  }
  elapsed = (diff.tv_sec)*MILLION + (diff.tv_nsec)/1000L;
  return elapsed;
}
#endif

void chefio_addReadBytes(size_t b) {
  const int i = chefio_global_stats.fileIdx;
  chefio_global_stats.readBytes[i] += b;
}

void chefio_addWriteBytes(size_t b) {
  const int i = chefio_global_stats.fileIdx;
  chefio_global_stats.writeBytes[i] += b;
}

void chefio_addReadTime(size_t t) {
  const int i = chefio_global_stats.fileIdx;
  chefio_global_stats.readTime[i] += t;
  chefio_global_stats.reads[i]++;
}

void chefio_addWriteTime(size_t t) {
  const int i = chefio_global_stats.fileIdx;
  chefio_global_stats.writeTime[i] += t;
  chefio_global_stats.writes[i]++;
}

void chefio_setfile(int f) {
  assert(f >= 0 && f < NUM_CHEF_FILES);
  chefio_global_stats.fileIdx = f;
}

void chefio_addOpenTime(size_t t) {
  const int i = chefio_global_stats.fileIdx;
  chefio_global_stats.openTime[i] += t;
  chefio_global_stats.opens[i]++;
}

void chefio_addCloseTime(size_t t) {
  const int i = chefio_global_stats.fileIdx;
  chefio_global_stats.closeTime[i] += t;
  chefio_global_stats.closes[i]++;
}

static void printMinMaxAvgSzt(const char* key, size_t v) {
  size_t min = PCU_Min_SizeT(v);
  size_t max = PCU_Max_SizeT(v);
  size_t tot = PCU_Add_SizeT(v);
  double avg = ((double)tot)/PCU_Comm_Peers();
  if(!PCU_Comm_Self())
    fprintf(stderr, "chefio_%s min max avg %" PRIu64 " %" PRIu64 " %f\n", key, min, max, avg);
}

static void printMinMaxAvgDbl(const char* key, double v) {
  double min = PCU_Min_Double(v);
  double max = PCU_Max_Double(v);
  double tot = PCU_Add_Double(v);
  double avg = tot/PCU_Comm_Peers();
  if(!PCU_Comm_Self())
    fprintf(stderr, "chefio_%s min max avg %f %f %f\n",
        key, min, max, avg);
}

static size_t chefio_getReadTime(int file) {
  return chefio_global_stats.readTime[file];
}

static size_t chefio_getWriteTime(int file) {
  return chefio_global_stats.writeTime[file];
}

static size_t chefio_getReadBytes(int file) {
  return chefio_global_stats.readBytes[file];
}

static size_t chefio_getWriteBytes(int file) {
  return chefio_global_stats.writeBytes[file];
}

static size_t chefio_getReads(int file) {
  return chefio_global_stats.reads[file];
}

static size_t chefio_getWrites(int file) {
  return chefio_global_stats.writes[file];
}

static size_t chefio_getOpens(int file) {
  return chefio_global_stats.opens[file];
}

static size_t chefio_getCloses(int file) {
  return chefio_global_stats.closes[file];
}

static size_t chefio_getOpenTime(int file) {
  return chefio_global_stats.openTime[file];
}

static size_t chefio_getCloseTime(int file) {
  return chefio_global_stats.closeTime[file];
}

void chefio_printStats() {
  if(!PCU_Comm_Self()) {
    const size_t us = 1000;
    chefioTime t0,t1;
    size_t elapsed;
    chefio_time(&t0);
    usleep(us);
    chefio_time(&t1);
    elapsed = chefio_time_diff(&t0,&t1);
    fprintf(stderr, "%" PRIu64 " us measured as %" PRIu64 " us\n", us, elapsed);
  }
  for(int chefFile=0; chefFile<NUM_CHEF_FILES; chefFile++) {
    if(!PCU_Comm_Self())
      fprintf(stderr, "chefio_filename %s\n", chefFile ? "restart" : "geombc");
    int reads = PCU_Max_Int((int)chefio_getReads(chefFile));
    if(reads) {
      printMinMaxAvgSzt("reads", chefio_getReads(chefFile));
      printMinMaxAvgSzt("readTime (us)", chefio_getReadTime(chefFile));
      printMinMaxAvgSzt("readBytes (B)", chefio_getReadBytes(chefFile));
      double bw = ((double)chefio_getReadBytes(chefFile))/chefio_getReadTime(chefFile);
      printMinMaxAvgDbl("readBandwidth (MB/s)", bw);
      /* B  * 10^6us *  1MB   = MB
       * -    ------   -----    --
       * us     1s     10^6B    s
       */
    }
    int writes = PCU_Max_Int((int)chefio_getWrites(chefFile));
    if(writes) {
      printMinMaxAvgSzt("writes", chefio_getWrites(chefFile));
      printMinMaxAvgSzt("writeTime (us)", chefio_getWriteTime(chefFile)); //FIXME
      printMinMaxAvgSzt("writeBytes (B)", chefio_getWriteBytes(chefFile));
      printMinMaxAvgDbl("writeBandwidth (MB/s)",
          ((double)chefio_getWriteBytes(chefFile))/chefio_getWriteTime(chefFile));
    }
    int opens = PCU_Max_Int((int)chefio_getOpens(chefFile));
    if(opens) {
      printMinMaxAvgSzt("opens", chefio_getOpens(chefFile));
      printMinMaxAvgSzt("openTime (us)", chefio_getOpenTime(chefFile));
    }
    int closes = PCU_Max_Int((int)chefio_getCloses(chefFile));
    if(closes) {
      printMinMaxAvgSzt("closes", chefio_getCloses(chefFile));
      printMinMaxAvgSzt("closeTime (us)", chefio_getCloseTime(chefFile));
    }
  }
}

void chefio_initStats() {
#ifdef __INTEL_COMPILER
  chefio_global_stats.cpus = chefio_getCyclesPerMicroSec();
#endif
  for(int i=0; i<NUM_CHEF_FILES; i++) {
    chefio_global_stats.readTime[i] = 0;
    chefio_global_stats.writeTime[i] = 0;
    chefio_global_stats.readBytes[i] = 0;
    chefio_global_stats.writeBytes[i] = 0;
    chefio_global_stats.reads[i] = 0;
    chefio_global_stats.writes[i] = 0;
    chefio_global_stats.openTime[i] = 0;
    chefio_global_stats.closeTime[i] = 0;
    chefio_global_stats.opens[i] = 0;
    chefio_global_stats.closes[i] = 0;
  }
}

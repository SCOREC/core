#include <stdio.h>
#include <pcu_util.h>
#include <PCU.h>
#include <PCUObj.h>
#include <phiotimer.h>
#include <lionPrint.h>
#include <sstream>

#include <time.h> /* clock_gettime */
#include <unistd.h> /* usleep */

#define BILLION 1000L*1000L*1000L
#define MILLION 1000L*1000L
struct phastaio_stats {
  size_t cpus;
  size_t readTime[NUM_PHASTAIO_MODES];
  size_t writeTime[NUM_PHASTAIO_MODES];
  size_t readBytes[NUM_PHASTAIO_MODES];
  size_t writeBytes[NUM_PHASTAIO_MODES];
  size_t reads[NUM_PHASTAIO_MODES];
  size_t writes[NUM_PHASTAIO_MODES];
  size_t openTime[NUM_PHASTAIO_MODES];
  size_t closeTime[NUM_PHASTAIO_MODES];
  size_t opens[NUM_PHASTAIO_MODES];
  size_t closes[NUM_PHASTAIO_MODES];
  int fileIdx;
};
static struct phastaio_stats phastaio_global_stats;

#if defined(HAVE_INTEL_RDTSC)
/* return the cycle count */
void phastaio_time(phastaioTime* t) {
  *t = _rdtsc(); //intel intrinsic
}
/* determine the reference clock frequency */
static size_t phastaio_getCyclesPerMicroSec(PCU_t h) {
  const size_t usec = 5*MILLION;
  size_t cpus, cycles;
  phastaioTime t0, t1;
  phastaio_time(&t0);
  /* Testing on Theta indicates that 5s is long enough
   * to get a stable value for the reference frequency.
   */
  usleep(usec);
  phastaio_time(&t1);
  cycles = t1 - t0;
  cpus = ((double)cycles)/(usec);
  if(!PCU_Comm_Self(h)) {
    std::stringstream ss;
    ss << "cycles " << cycles << " us " << usec
       << " cycles per micro second " << cpus << "\n";
    std::string s = ss.str();
    lion_eprint(1,"%s",s.c_str());
  }
  return cpus;
}
/*return elapsed time in micro seconds*/
size_t phastaio_time_diff(phastaioTime* start, phastaioTime* end) {
  size_t cycles = *end - *start;
  size_t us = ((double)cycles)/phastaio_global_stats.cpus;
  return us;
}
#elif defined(USE_PCU_TIME)
void phastaio_time(phastaioTime* t) {
  *t = pcu::Time();
}
/*return elapsed time in micro seconds*/
size_t phastaio_time_diff(phastaioTime* start, phastaioTime* end) {
  size_t elapsed = static_cast<size_t>((*end-*start)*MILLION);
  return elapsed;
}
#elif defined(HAVE_CLOCK_GETTIME)
void phastaio_time(phastaioTime* t) {
  int err;
  err = clock_gettime(CLOCK_MONOTONIC,t);
  PCU_ALWAYS_ASSERT(!err);
}
/*return elapsed time in micro seconds*/
size_t phastaio_time_diff(phastaioTime* start, phastaioTime* end) {
  PCU_ALWAYS_ASSERT(sizeof(size_t)==8);
  size_t elapsed = 0;
  phastaioTime diff;
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

void phastaio_addReadBytes(size_t b) {
  const int i = phastaio_global_stats.fileIdx;
  phastaio_global_stats.readBytes[i] += b;
}

void phastaio_addWriteBytes(size_t b) {
  const int i = phastaio_global_stats.fileIdx;
  phastaio_global_stats.writeBytes[i] += b;
}

void phastaio_addReadTime(size_t t) {
  const int i = phastaio_global_stats.fileIdx;
  phastaio_global_stats.readTime[i] += t;
  phastaio_global_stats.reads[i]++;
}

void phastaio_addWriteTime(size_t t) {
  const int i = phastaio_global_stats.fileIdx;
  phastaio_global_stats.writeTime[i] += t;
  phastaio_global_stats.writes[i]++;
}

void phastaio_setfile(int f) {
  char msg[64]; sprintf(msg, "f %d", f);
  PCU_ALWAYS_ASSERT_VERBOSE(f >= 0 && f < NUM_PHASTAIO_MODES, msg);
  phastaio_global_stats.fileIdx = f;
}

void phastaio_addOpenTime(size_t t) {
  const int i = phastaio_global_stats.fileIdx;
  phastaio_global_stats.openTime[i] += t;
  phastaio_global_stats.opens[i]++;
}

void phastaio_addCloseTime(size_t t) {
  const int i = phastaio_global_stats.fileIdx;
  phastaio_global_stats.closeTime[i] += t;
  phastaio_global_stats.closes[i]++;
}

static const char* getFileName() {
  const char* names[NUM_PHASTAIO_MODES] = {
    "geombc_read",
    "geombc_write",
    "restart_read",
    "restart_write"};
  return names[phastaio_global_stats.fileIdx];
}

static void printMinMaxAvgSzt(const char* key, size_t v, PCU_t h) {
  size_t min = PCU_Min_SizeT(h, v);
  size_t max = PCU_Max_SizeT(h, v);
  size_t tot = PCU_Add_SizeT(h, v);
  double avg = ((double)tot)/PCU_Comm_Peers(h);
  if(!PCU_Comm_Self(h)) {
    std::stringstream ss;
    ss << getFileName() << "_" << key << "min max avg"
       << min << " " << max << " " << avg << "\n";
    std::string s = ss.str();
    lion_eprint(1,"%s",s.c_str());
  }
}

static void printMinMaxAvgDbl(const char* key, double v, PCU_t h) {
  double min = PCU_Min_Double(h, v);
  double max = PCU_Max_Double(h, v);
  double tot = PCU_Add_Double(h, v);
  double avg = tot/PCU_Comm_Peers(h);
  if(!PCU_Comm_Self(h))
    lion_eprint(1, "%s_%s min max avg %f %f %f\n",
        getFileName(), key, min, max, avg);
}

static size_t phastaio_getReadTime() {
  const int i = phastaio_global_stats.fileIdx;
  return phastaio_global_stats.readTime[i];
}

static size_t phastaio_getWriteTime() {
  const int i = phastaio_global_stats.fileIdx;
  return phastaio_global_stats.writeTime[i];
}

static size_t phastaio_getReadBytes() {
  const int i = phastaio_global_stats.fileIdx;
  return phastaio_global_stats.readBytes[i];
}

static size_t phastaio_getWriteBytes() {
  const int i = phastaio_global_stats.fileIdx;
  return phastaio_global_stats.writeBytes[i];
}

static size_t phastaio_getReads() {
  const int i = phastaio_global_stats.fileIdx;
  return phastaio_global_stats.reads[i];
}

static size_t phastaio_getWrites() {
  const int i = phastaio_global_stats.fileIdx;
  return phastaio_global_stats.writes[i];
}

static size_t phastaio_getOpens() {
  const int i = phastaio_global_stats.fileIdx;
  return phastaio_global_stats.opens[i];
}

static size_t phastaio_getCloses() {
  const int i = phastaio_global_stats.fileIdx;
  return phastaio_global_stats.closes[i];
}

static size_t phastaio_getOpenTime() {
  const int i = phastaio_global_stats.fileIdx;
  return phastaio_global_stats.openTime[i];
}

static size_t phastaio_getCloseTime() {
  const int i = phastaio_global_stats.fileIdx;
  return phastaio_global_stats.closeTime[i];
}

void phastaio_printStats(PCU_t h) {
  if(!PCU_Comm_Self(h)) {
    const size_t us = 1000;
    phastaioTime t0,t1;
    size_t elapsed;
    phastaio_time(&t0);
    usleep(us);
    phastaio_time(&t1);
    elapsed = phastaio_time_diff(&t0,&t1);
    std::stringstream ss;
    ss << us << " us measured as " << elapsed << "us\n";
    std::string s = ss.str();
    lion_eprint(1,"%s",s.c_str());
  }
  for(int chefFile=0; chefFile<NUM_PHASTAIO_MODES; chefFile++) {
    size_t totalus = 0;
    size_t totalbytes = 0;
    phastaio_setfile(chefFile);
    if(!PCU_Comm_Self(h))
      lion_eprint(1, "phastaio_filename %s\n", getFileName());
    int reads = PCU_Max_Int(h, (int)phastaio_getReads());
    if(reads) {
      totalus += phastaio_getReadTime();
      totalbytes += phastaio_getReadBytes();
      printMinMaxAvgSzt("reads", phastaio_getReads(), h);
      printMinMaxAvgSzt("readTime (us)", phastaio_getReadTime(), h);
      printMinMaxAvgSzt("readBytes (B)", phastaio_getReadBytes(), h);
      double bw = ((double)phastaio_getReadBytes())/phastaio_getReadTime();
      printMinMaxAvgDbl("readBandwidth (MB/s)", bw, h);
      /* B  * 10^6us *  1MB   = MB
       * -    ------   -----    --
       * us     1s     10^6B    s
       */
    }
    int writes = PCU_Max_Int(h, (int)phastaio_getWrites());
    if(writes) {
      totalus += phastaio_getWriteTime();
      totalbytes += phastaio_getWriteBytes();
      printMinMaxAvgSzt("writes", phastaio_getWrites(), h);
      printMinMaxAvgSzt("writeTime (us)", phastaio_getWriteTime(), h);
      printMinMaxAvgSzt("writeBytes (B)", phastaio_getWriteBytes(), h);
      printMinMaxAvgDbl("writeBandwidth (MB/s)",
          ((double)phastaio_getWriteBytes())/phastaio_getWriteTime(), h);
    }
    int opens = PCU_Max_Int(h, (int)phastaio_getOpens());
    if(opens) {
      totalus += phastaio_getOpenTime();
      printMinMaxAvgSzt("opens", phastaio_getOpens(), h);
      printMinMaxAvgSzt("openTime (us)", phastaio_getOpenTime(), h);
    }
    int closes = PCU_Max_Int(h, (int)phastaio_getCloses());
    if(closes) {
      totalus += phastaio_getCloseTime();
      printMinMaxAvgSzt("closes", phastaio_getCloses(), h);
      printMinMaxAvgSzt("closeTime (us)", phastaio_getCloseTime(), h);
    }
    if(totalbytes) {
      printMinMaxAvgSzt("totalTime (us)", totalus, h);
      printMinMaxAvgSzt("totalBytes (B)", totalbytes, h);
      printMinMaxAvgDbl("effectiveBandwidth (MB/s)",
          ((double)totalbytes)/totalus, h);
    }
  }
}

void phastaio_initStats(PCU_t h) {
  (void) h;
#ifdef __INTEL_COMPILER
  phastaio_global_stats.cpus = phastaio_getCyclesPerMicroSec(h);
#endif
  for(int i=0; i<NUM_PHASTAIO_MODES; i++) {
    phastaio_global_stats.readTime[i] = 0;
    phastaio_global_stats.writeTime[i] = 0;
    phastaio_global_stats.readBytes[i] = 0;
    phastaio_global_stats.writeBytes[i] = 0;
    phastaio_global_stats.reads[i] = 0;
    phastaio_global_stats.writes[i] = 0;
    phastaio_global_stats.openTime[i] = 0;
    phastaio_global_stats.closeTime[i] = 0;
    phastaio_global_stats.opens[i] = 0;
    phastaio_global_stats.closes[i] = 0;
  }
}

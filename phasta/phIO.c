#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <pcu_io.h>
#include <phIO.h>
#include <PCU.h>
#include <inttypes.h> /* PRIu64 */
#include <time.h> /* clock_gettime */
#include <unistd.h> /* usleep */

#define PH_LINE 1024
#define MAGIC 362436
#define FIELD_PARAMS 3

#define BILLION 1000L*1000L*1000L
#define MILLION 1000L*1000L

#ifdef __INTEL_COMPILER
typedef size_t chefioTime;
size_t chefio_global_cpus;
static size_t chefio_time_diff(chefioTime* start, chefioTime* end);
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
static size_t chefio_time_diff(chefioTime* start, chefioTime* end) {
  size_t cycles = *end - *start;
  size_t us = ((double)cycles)/chefio_global_cpus;
  return us;
}
#else
typedef struct timespec chefioTime;
static void chefio_time(chefioTime* t) {
  int err;
  err = clock_gettime(CLOCK_MONOTONIC,t);
  assert(!err);
}
/*return elapsed time in micro seconds*/
static size_t chefio_time_diff(chefioTime* start, chefioTime* end) {
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

struct chefio_stats {
  size_t readTime;
  size_t writeTime;
  size_t readBytes;
  size_t writeBytes;
  size_t reads;
  size_t writes;
};
static struct chefio_stats chefio_global_stats;

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

static size_t chefio_getReadTime() {
  return chefio_global_stats.readTime;
}

static size_t chefio_getWriteTime() {
  return chefio_global_stats.writeTime;
}

static size_t chefio_getReadBytes() {
  return chefio_global_stats.readBytes;
}

static size_t chefio_getWriteBytes() {
  return chefio_global_stats.writeBytes;
}

static size_t chefio_getReads() {
  return chefio_global_stats.reads;
}

static size_t chefio_getWrites() {
  return chefio_global_stats.writes;
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
  int reads = PCU_Max_Int((int)chefio_getReads());
  if(reads) {
    printMinMaxAvgSzt("reads", chefio_getReads());
    printMinMaxAvgSzt("readTime (us)", chefio_getReadTime());
    printMinMaxAvgSzt("readBytes (B)", chefio_getReadBytes());
    double bw = ((double)chefio_getReadBytes())/chefio_getReadTime();
    if( chefio_getReadTime() == 0 ) {
      fprintf(stderr, "%d ZERO read time reads %" PRIu64 " readTime %" PRIu64 " readBytes %" PRIu64 " bw %f\n",
          PCU_Comm_Self(), chefio_getReads(), chefio_getReadTime(), chefio_getReadBytes(), bw);
    }
    printMinMaxAvgDbl("readBandwidth (MB/s)", bw);
    /* B  * 10^6us *  1MB   = MB
     * -    ------   -----    --
     * us     1s     10^6B    s
     */
  }
  int writes = PCU_Max_Int((int)chefio_getWrites());
  if(writes) {
    printMinMaxAvgSzt("writes", chefio_getWrites());
    printMinMaxAvgSzt("writeTime (us)", chefio_getWriteTime());
    printMinMaxAvgSzt("writeBytes (B)", chefio_getWriteBytes());
    printMinMaxAvgDbl("writeBandwidth (MB/s)",
        ((double)chefio_getWriteBytes())/chefio_getWriteTime());
  }
}

void chefio_initStats() {
#ifdef __INTEL_COMPILER
  chefio_global_cpus = chefio_getCyclesPerMicroSec();
#endif
  chefio_global_stats.readTime = 0;
  chefio_global_stats.writeTime = 0;
  chefio_global_stats.readBytes = 0;
  chefio_global_stats.writeBytes = 0;
  chefio_global_stats.reads = 0;
  chefio_global_stats.writes = 0;
}

enum {
NODES_PARAM,
VARS_PARAM,
STEP_PARAM
};

static const char* magic_name = "byteorder magic number";

void ph_write_header(FILE* f, const char* name, size_t bytes,
    int nparam, int* params)
{
  int i;
  fprintf(f,"%s : < %lu > ", name, (long)bytes);
  for (i = 0; i < nparam; ++i)
    fprintf(f, "%d ", params[i]);
  fprintf(f, "\n");
}

static void skip_leading_spaces(char** s)
{
  while (**s == ' ') ++(*s);
}

static void cut_trailing_spaces(char* s)
{
  char* e = s + strlen(s);
  for (--e; e >= s; --e)
    if (*e != ' ')
      break;
  ++e;
  *e = '\0';
}

static void parse_header(char* header, char** name, long* bytes,
    int nparam, int* params)
{
  char* saveptr = NULL;
  int i = 0;
  assert(header != NULL);
  header = strtok_r(header, ":", &saveptr);
  if (name) {
    *name = header;
    skip_leading_spaces(name);
    cut_trailing_spaces(*name);
  }
  strtok_r(NULL, "<", &saveptr);
  header = strtok_r(NULL, ">", &saveptr);
  if (bytes)
    sscanf(header, "%ld", bytes);
  if (params) {
    while( (header = strtok_r(NULL, " \n", &saveptr)) )
      sscanf(header, "%d", &params[i++]);
    while( i < nparam )
      params[i++] = 0;
  }
}

static int find_header(FILE* f, const char* name, char* found, char header[PH_LINE])
{
  char* hname;
  long bytes;
  char tmp[PH_LINE];
  while (fgets(header, PH_LINE, f)) {
    if ((header[0] == '#') || (header[0] == '\n'))
      continue;
    strncpy(tmp, header, PH_LINE-1);
    tmp[PH_LINE-1] = '\0';
    parse_header(tmp, &hname, &bytes, 0, NULL);
    if (!strncmp(name, hname, strlen(name))) {
      strncpy(found, hname, strlen(hname));
      found[strlen(hname)] = '\0';
      return 1;
    }
    fseek(f, bytes, SEEK_CUR);
  }
  if (!PCU_Comm_Self() && strlen(name) > 0)
    fprintf(stderr,"warning: phIO could not find \"%s\"\n",name);
  return 0;
}

static void write_magic_number(FILE* f)
{
  int why = 1;
  int magic = MAGIC;
  ph_write_header(f, magic_name, sizeof(int) + 1, 1, &why);
  fwrite(&magic, sizeof(int), 1, f);
  fprintf(f,"\n");
}

static int seek_after_header(FILE* f, const char* name)
{
  char dummy[PH_LINE];
  char found[PH_LINE];
  return find_header(f, name, found, dummy);
}

static void my_fread(void* p, size_t size, size_t nmemb, FILE* f)
{
  chefioTime t0,t1;
  chefio_time(&t0);
  size_t r = fread(p, size, nmemb, f);
  chefio_time(&t1);
  const size_t time = chefio_time_diff(&t0,&t1);
  const size_t bytes = nmemb*size;
  chefio_global_stats.readTime += time;
  chefio_global_stats.readBytes += bytes;
  chefio_global_stats.reads++;
  if( !time ) {
    fprintf(stderr, "%s %d ZERO read time reads %" PRIu64 " time %" PRIu64 " bytes %" PRIu64 "\n",
        __func__, PCU_Comm_Self(), chefio_getReads(), time, bytes);
  }
  assert(r == nmemb);
}

static int read_magic_number(FILE* f)
{
  int magic;
  if (!seek_after_header(f, magic_name)) {
    if (!PCU_Comm_Self())
      fprintf(stderr,"warning: not swapping bytes\n");
    rewind(f);
    return 0;
  }
  my_fread(&magic, sizeof(int), 1, f);
  return magic != MAGIC;
}

void ph_write_preamble(FILE* f)
{
  fprintf(f, "# PHASTA Input File Version 2.0\n");
  fprintf(f, "# Byte Order Magic Number : 362436 \n");
  fprintf(f, "# Output generated by libph version: yes\n");
  write_magic_number(f);
}

void ph_write_doubles(FILE* f, const char* name, double* data,
    size_t n, int nparam, int* params)
{
  chefioTime t0,t1;
  ph_write_header(f, name, n * sizeof(double) + 1, nparam, params);
  chefio_time(&t0);
  fwrite(data, sizeof(double), n, f);
  chefio_time(&t1);
  fprintf(f, "\n");
  chefio_global_stats.writeTime += chefio_time_diff(&t0,&t1);
  chefio_global_stats.writeBytes += n*sizeof(double);
  chefio_global_stats.writes++;
}

void ph_write_ints(FILE* f, const char* name, int* data,
    size_t n, int nparam, int* params)
{
  chefioTime t0,t1;
  ph_write_header(f, name, n * sizeof(int) + 1, nparam, params);
  chefio_time(&t0);
  fwrite(data, sizeof(int), n, f);
  chefio_time(&t1);
  fprintf(f, "\n");
  chefio_global_stats.writeTime += chefio_time_diff(&t0,&t1);
  chefio_global_stats.writeBytes += n*sizeof(int);
  chefio_global_stats.writes++;
}

static void parse_params(char* header, long* bytes,
    int* nodes, int* vars, int* step)
{
  int params[FIELD_PARAMS];
  parse_header(header, NULL, bytes, FIELD_PARAMS, params);
  *nodes = params[NODES_PARAM];
  *vars = params[VARS_PARAM];
  *step = params[STEP_PARAM];
}

int ph_should_swap(FILE* f) {
  return read_magic_number(f);
}

int ph_read_field(FILE* f, const char* field, int swap,
    double** data, int* nodes, int* vars, int* step, char* hname)
{
  long bytes, n;
  char header[PH_LINE];
  int ok;
  ok = find_header(f, field, hname, header);
  if(!ok) /* not found */
    return 0;
  parse_params(header, &bytes, nodes, vars, step);
  if(!bytes) /* empty data block */
    return 1;
  assert(((bytes - 1) % sizeof(double)) == 0);
  n = (bytes - 1) / sizeof(double);
  assert((int)n == (*nodes) * (*vars));
  *data = malloc(bytes);
  my_fread(*data, sizeof(double), n, f);
  if (swap)
    pcu_swap_doubles(*data, n);
  return 2;
}

void ph_write_field(FILE* f, const char* field, double* data,
    int nodes, int vars, int step)
{
  int params[FIELD_PARAMS];
  params[NODES_PARAM] = nodes;
  params[VARS_PARAM] = vars;
  params[STEP_PARAM] = step;
  ph_write_doubles(f, field, data, nodes * vars, FIELD_PARAMS, params);
}

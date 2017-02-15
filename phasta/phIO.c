#include <assert.h>
#include <string.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <pcu_io.h>
#include <phIO.h>
#include <PCU.h>

#define PH_LINE 1024
#define MAGIC 362436
#define FIELD_PARAMS 3

struct chefio_stats {
  double readTime;
  double writeTime;
  size_t readBytes;
  size_t writeBytes;
  size_t reads;
  size_t writes;
};
struct chefio_stats chefio_global_stats;

void printMinMaxAvgSzt(const char* key, size_t v) {
  int val = (int)v;
  int min = PCU_Min_Int(val);
  int max = PCU_Max_Int(val);
  long tot = PCU_Add_Long((long)val);
  double avg = tot/(double)PCU_Comm_Peers();
  if(!PCU_Comm_Self())
    fprintf(stderr, "chefio_%s min max avg %d %d %f\n",
        key, min, max, avg);
}

void printMinMaxAvgDbl(const char* key, double v) {
  double min = PCU_Min_Double(v);
  double max = PCU_Max_Double(v);
  double tot = PCU_Add_Double(v);
  double avg = tot/PCU_Comm_Peers();
  if(!PCU_Comm_Self())
    fprintf(stderr, "chefio_%s min max avg %f %f %f\n",
        key, min, max, avg);
}

double chefio_getReadTime() {
  return chefio_global_stats.readTime;
}

double chefio_getWriteTime() {
  return chefio_global_stats.writeTime;
}

size_t chefio_getReadBytes() {
  return chefio_global_stats.readBytes;
}

size_t chefio_getWriteBytes() {
  return chefio_global_stats.writeBytes;
}

size_t chefio_getReads() {
  return chefio_global_stats.reads;
}

size_t chefio_getWrites() {
  return chefio_global_stats.writes;
}

void chefio_printStats() {
  const int mebi=1024*1024;
  int reads = PCU_Max_Int((int)chefio_getReads());
  if(reads) {
    printMinMaxAvgDbl("readTime (s)",chefio_getReadTime());
    printMinMaxAvgSzt("readBytes (B)", chefio_getReadBytes());
    printMinMaxAvgDbl("readBandwidth (MiB/s)",
        (chefio_getReadBytes()/chefio_getReadTime())/mebi);
  }
  int writes = PCU_Max_Int((int)chefio_getWrites());
  if(writes) {
    printMinMaxAvgDbl("writeTime (s)", chefio_getWriteTime());
    printMinMaxAvgSzt("writeBytes (B)", chefio_getWriteBytes());
    printMinMaxAvgDbl("writeBandwidth (MiB/s)",
        (chefio_getWriteBytes()/chefio_getWriteTime())/mebi);
  }
}

void chefio_initStats() {
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
  double t0 = PCU_Time();
  int i;
  fprintf(f,"%s : < %lu > ", name, (long)bytes);
  for (i = 0; i < nparam; ++i)
    fprintf(f, "%d ", params[i]);
  fprintf(f, "\n");
  chefio_global_stats.writeTime += PCU_Time()-t0;
  chefio_global_stats.writeBytes += nparam*sizeof(int);
  chefio_global_stats.writes++;
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
  double t0 = PCU_Time();
  size_t r = fread(p, size, nmemb, f);
  assert(r == nmemb);
  chefio_global_stats.readTime += PCU_Time()-t0;
  chefio_global_stats.readBytes += nmemb*size;
  chefio_global_stats.reads++;
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
  double t0 = PCU_Time();
  ph_write_header(f, name, n * sizeof(double) + 1, nparam, params);
  fwrite(data, sizeof(double), n, f);
  fprintf(f, "\n");
  chefio_global_stats.writeTime += PCU_Time()-t0;
  chefio_global_stats.writeBytes += n*sizeof(double);
  chefio_global_stats.writes++;
}

void ph_write_ints(FILE* f, const char* name, int* data,
    size_t n, int nparam, int* params)
{
  double t0 = PCU_Time();
  ph_write_header(f, name, n * sizeof(int) + 1, nparam, params);
  fwrite(data, sizeof(int), n, f);
  fprintf(f, "\n");
  chefio_global_stats.writeTime += PCU_Time()-t0;
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

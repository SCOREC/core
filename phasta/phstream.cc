#include <stdlib.h>
#include <string>
#include <assert.h>
#include "phstream.h"
#include <mpi.h>

#ifndef PHSTREAM_TIMERS_ON
#define PHSTREAM_TIMERS_ON 0
#endif

namespace {
  inline double getTime() {
    return MPI_Wtime();
  }
  inline bool isRankZero() {
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return !rank;
  }
  inline void printTime(const char* key, double t) {
#if PHSTREAM_TIMERS_ON==1
    if( isRankZero() )
      fprintf(stderr, "%s %f seconds\n", key, t);
#endif
  }
}

extern "C" {
  struct RStream{
    char* restart;
    size_t rSz;
  };
}

extern "C" {
  struct GRStream{
    char *geom, *restart;
    size_t gSz, rSz;
  };
}

RStream* makeRStream() {
  const double t0 = getTime();
  RStream* rs = (RStream*) malloc(sizeof(RStream));
  rs->restart = NULL;
  rs->rSz = 0;
  printTime(__func__, getTime()-t0);
  return rs;
}

FILE* openRStreamRead(RStream* rs) {
  const double t0 = getTime();
  FILE* f = fmemopen(rs->restart, rs->rSz, "r");
  printTime(__func__, getTime()-t0);
  return f;
}

FILE* openRStreamWrite(RStream* rs) {
  const double t0 = getTime();
  FILE* f = open_memstream(&(rs->restart), &(rs->rSz));
  printTime(__func__, getTime()-t0);
  return f;
}

void clearRStream(RStream* rs) {
  const double t0 = getTime();
  if(rs->restart) {
    free(rs->restart);
    rs->restart = NULL;
    rs->rSz = 0;
  }
  printTime(__func__, getTime()-t0);
}

void destroyRStream(RStream* rs) {
  const double t0 = getTime();
  clearRStream(rs);
  free(rs);
  printTime(__func__, getTime()-t0);
}

void attachRStream(GRStream* grs, RStream* rs) {
  const double t0 = getTime();
  rs->restart = grs->restart;
  rs->rSz = grs->rSz;
  grs->restart = NULL;
  grs->rSz = 0;
  printTime(__func__, getTime()-t0);
}


GRStream* makeGRStream() {
  const double t0 = getTime();
  GRStream* grs = (GRStream*) malloc(sizeof(GRStream));
  grs->geom = NULL;
  grs->gSz = 0;
  grs->restart = NULL;
  grs->rSz = 0;
  printTime(__func__, getTime()-t0);
  return grs;
}

void whichStream(const char* name, bool& isR, bool& isG) {
  const double t0 = getTime();
  std::string fname(name);
  std::string restartStr("restart");
  std::string geombcStr("geombc");
  isR = (fname.find(restartStr) != std::string::npos);
  isG = (fname.find(geombcStr)  != std::string::npos);
  assert(isR != isG);
  printTime(__func__, getTime()-t0);
}

void writeUnknown(const char* fname) {
  fprintf(stderr,
      "ERROR %s type of stream %s is unknown... exiting\n",
      __func__, fname);
}

FILE* openGRStreamRead(GRStream* grs, const char* named) {
  const double t0 = getTime();
  bool isR, isG;
  whichStream(named, isR, isG);
  FILE* f = NULL;
  if( isR && !isG )
    f = fmemopen(grs->restart, grs->rSz, "r");
  else if( isG && !isR )
    f = fmemopen(grs->geom, grs->gSz, "r");
  else {
    writeUnknown(named);
    exit(1);
  }
  printTime(__func__, getTime()-t0);
  return f;
}

FILE* openGRStreamWrite(GRStream* grs, const char* named) {
  const double t0 = getTime();
  bool isR, isG;
  whichStream(named, isR, isG);
  FILE* f = NULL;
  if( isR && !isG )
    f = open_memstream(&(grs->restart), &(grs->rSz));
  else if( isG && !isR )
    f = open_memstream(&(grs->geom), &(grs->gSz));
  else {
    writeUnknown(named);
    exit(1);
  }
  printTime(__func__, getTime()-t0);
  return f;
}

void clearGRStream(GRStream* grs) {
  const double t0 = getTime();
  if(grs->geom) {
    free(grs->geom);
    grs->geom = NULL;
    grs->gSz = 0;
  }
  if(grs->restart) {
    free(grs->restart);
    grs->restart = NULL;
    grs->rSz = 0;
  }
  printTime(__func__, getTime()-t0);
}

void destroyGRStream(GRStream* grs) {
  const double t0 = getTime();
  clearGRStream(grs);
  free(grs);
  printTime(__func__, getTime()-t0);
}


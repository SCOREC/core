#include <stdlib.h>
#include <string>
#include <pcu_util.h>
#include <lionPrint.h>
#include "phstream.h"
//#include <PCU.h>

#ifndef PHSTREAM_TIMERS_ON
#define PHSTREAM_TIMERS_ON 0
#endif

namespace {
  inline double getTime() {
    return pcu::Time();
  }
#if PHSTREAM_TIMERS_ON==1
  inline bool isRankZero(pcu::PCU *pcu_obj) {
    int rank = 0;
    rank = pcu_obj->Self();
    return !rank;
  }
#endif
  inline void printTime(const char* key, double t, pcu::PCU *pcu_obj) {
    (void) key;
    (void) t;
    (void) pcu_obj;
#if PHSTREAM_TIMERS_ON==1
    if( isRankZero(pcu_obj) )
      lion_eprint(1, "%s %f seconds\n", key, t);
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

RStream* makeRStream(pcu::PCU *PCUObj) {
  const double t0 = getTime();
  RStream* rs = (RStream*) malloc(sizeof(RStream));
  rs->restart = NULL;
  rs->rSz = 0;
  printTime(__func__, getTime()-t0, PCUObj);
  return rs;
}

#ifdef __APPLE__
FILE* openRStreamRead(RStream*, pcu::PCU*) {
  return NULL;
}
#else
FILE* openRStreamRead(RStream* rs, pcu::PCU *PCUObj) {
  const double t0 = getTime();
  FILE* f = fmemopen(rs->restart, rs->rSz, "r");
  printTime(__func__, getTime()-t0, PCUObj);
  return f;
}
#endif

#ifdef __APPLE__
FILE* openRStreamWrite(RStream*, pcu::PCU*) {
  return NULL;
}
#else
FILE* openRStreamWrite(RStream* rs, pcu::PCU *PCUObj) {
  const double t0 = getTime();
  FILE* f = open_memstream(&(rs->restart), &(rs->rSz));
  printTime(__func__, getTime()-t0, PCUObj);
  return f;
}
#endif

void clearRStream(RStream* rs, pcu::PCU *PCUObj) {
  const double t0 = getTime();
  if(rs->restart) {
    free(rs->restart);
    rs->restart = NULL;
    rs->rSz = 0;
  }
  printTime(__func__, getTime()-t0, PCUObj);
}

void destroyRStream(RStream* rs, pcu::PCU *PCUObj) {
  const double t0 = getTime();
  clearRStream(rs, PCUObj);
  free(rs);
  printTime(__func__, getTime()-t0, PCUObj);
}

void attachRStream(GRStream* grs, RStream* rs, pcu::PCU *PCUObj) {
  const double t0 = getTime();
  rs->restart = grs->restart;
  rs->rSz = grs->rSz;
  grs->restart = NULL;
  grs->rSz = 0;
  printTime(__func__, getTime()-t0, PCUObj);
}


GRStream* makeGRStream(pcu::PCU *PCUObj) {
  const double t0 = getTime();
  GRStream* grs = (GRStream*) malloc(sizeof(GRStream));
  grs->geom = NULL;
  grs->gSz = 0;
  grs->restart = NULL;
  grs->rSz = 0;
  printTime(__func__, getTime()-t0, PCUObj);
  return grs;
}

void whichStream(const char* name, bool& isR, bool& isG, pcu::PCU *PCUObj) {
  const double t0 = getTime();
  std::string fname(name);
  std::string restartStr("restart");
  std::string geombcStr("geombc");
  isR = (fname.find(restartStr) != std::string::npos);
  isG = (fname.find(geombcStr)  != std::string::npos);
  PCU_ALWAYS_ASSERT(isR != isG);
  printTime(__func__, getTime()-t0, PCUObj);
}

void writeUnknown(const char* fname) {
  lion_eprint(1,
      "ERROR %s type of stream %s is unknown... exiting\n",
      __func__, fname);
}

#ifdef __APPLE__
FILE* openGRStreamRead(GRStream*, const char*, pcu::PCU*) {
  return NULL;
}
#else
FILE* openGRStreamRead(GRStream* grs, const char* named, pcu::PCU *PCUObj) {
  const double t0 = getTime();
  bool isR, isG;
  whichStream(named, isR, isG, PCUObj);
  FILE* f = NULL;
  if( isR && !isG )
    f = fmemopen(grs->restart, grs->rSz, "r");
  else if( isG && !isR )
    f = fmemopen(grs->geom, grs->gSz, "r");
  else {
    writeUnknown(named);
    exit(1);
  }
  printTime(__func__, getTime()-t0, PCUObj);
  return f;
}
#endif

#ifdef __APPLE__
FILE* openGRStreamWrite(GRStream*, const char*, pcu::PCU*) {
  return NULL;
}
#else
FILE* openGRStreamWrite(GRStream* grs, const char* named, pcu::PCU *PCUObj) {
  const double t0 = getTime();
  bool isR, isG;
  whichStream(named, isR, isG, PCUObj);
  FILE* f = NULL;
  if( isR && !isG )
    f = open_memstream(&(grs->restart), &(grs->rSz));
  else if( isG && !isR )
    f = open_memstream(&(grs->geom), &(grs->gSz));
  else {
    writeUnknown(named);
    exit(1);
  }
  printTime(__func__, getTime()-t0, PCUObj);
  return f;
}
#endif

void clearGRStream(GRStream* grs, pcu::PCU *PCUObj) {
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
  printTime(__func__, getTime()-t0, PCUObj);
}

void destroyGRStream(GRStream* grs, pcu::PCU *PCUObj) {
  const double t0 = getTime();
  clearGRStream(grs, PCUObj);
  free(grs);
  printTime(__func__, getTime()-t0, PCUObj);
}


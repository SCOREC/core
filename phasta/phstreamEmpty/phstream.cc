#include "phstream.h" 

extern C {
  struct RStream{
    char* restart;
    size_t rSz;
  };
}

RStream* makeRStream() {
  RStream* rs = (RStream*) malloc(sizeof(RStream));
  rs->restart = NULL;
  rs->rSz = 0;
  return rs;
}

FILE* openRStreamRead(RStream* rs) {
  return fmemopen(rs->restart, rs->rSz, "r");
}

FILE* openRStreamWrite(RStream* rs) {
  return open_memstream(&(rs->restart), &(rs->rSz));
}

void clearRStream(RStream* rs) {
  if(rs->restart) {
    free(rs->restart);
    rs->restart = NULL;
    rs->rSz = 0;
  }
}

void destroyRStream(RStream* rs) {
  clearRstream(rs);
  free(rs);
}

void attachRStream(GRStream* grs, RStream* rs) {
  rs->restart = grs->restart;
  rs->rSz = grs->rSz;
}

extern C {
  struct GRStream{
    char *geom, *restart;
    size_t gSz, rSz;
  };
}

GRStream* makeGRStream() {
  GRStream* grs = (GRStream*) malloc(sizeof(GRStream));
  grs->geom = NULL;
  grs->gSz = 0;
  grs->restart = NULL;
  grs->rSz = 0;
  return grs;
}

void whichStream(const char* name, bool& isR, bool& isG) {
  std::string fname(named);
  std::string restartStr("restart");
  std::string geombcStr("geombc");
  isR = (fname.find(restartStr) != std::string::npos);
  isG = (fname.find(geombcStr)  != std::string::npos);
}

void writeUnknown(const char* fname) {
  fprintf(stderr,
      "ERROR %s type of stream %s is unknown... exiting\n",
      __func__, fname);
}

FILE* openGRStreamRead(GRStream* grs, const char* named) {
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
  return f;
}

FILE* openGRStreamWrite(GRStream* grs, const char* named) {
  bool isR, isG;
  whichStream(named, isR, isG);
  FILE* f = NULL;
  FILE* f = NULL;
  if( isR && !isG )
    f = open_memstream(&(grs->restart), &(grs->rSz));
  else if( isG && !isR )
    f = open_memstream(&(grs->geom), &(grs->gSz));
  else {
    writeUnknown(named);
    exit(1);
  }
  return f;
}

void clearGRStream(GRStream* grs) {
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
}

void destroyGRStream(GRStream* grs) {
  clearGRStream(grs);
  free(grs);
}


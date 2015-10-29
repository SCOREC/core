#include "phstream.h" 
#include <stdlib.h>
void fail(const char* f) {
  fprintf(stderr, 
    "ERROR: function %s is disabled - compile with chefPhasta\n", f);
  exit(EXIT_FAILURE);
}
FILE* openRStreamRead(RStream*) {
  fail(__func__);
  return NULL;
}
FILE* openGRStreamWrite(GRStream*, const char*) {
  fail(__func__);
  return NULL;
}


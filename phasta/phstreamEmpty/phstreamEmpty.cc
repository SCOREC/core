#include "phstream.h" 
void* fail(const char* f) {
  fprintf(stderr, 
    "ERROR: function %s is disabled - compile with chefPhasta\n", f);
  return NULL;
}
FILE* openRStreamRead(RStream*) {
  return fail(__func__);
}
FILE* openGRStreamWrite(GRStream*, const char*) {
  return fail(__func__);
}


#include "parmaIO.h"
#include "PCU.h"
#include <stdio.h>

void ParMA_IO::debugPrint(const char* msg, bool dbgOverride){
  if (dbgOverride) {
    printf("ParMA_DEBUG %s", msg);
    fflush(stdout);
  }    
}

void ParMA_IO::tracePrint(const char* fnName, bool traceOverride){
  if (traceOverride) 
    printf("[%d] Entering %s\n", PCU_Comm_Self(), fnName);
}

void ParMA_IO::statusPrint(const char* msg) {
  if ( 0 == PCU_Comm_Self() ) {
    printf("STATUS %s", msg);
    fflush(stdout);
  }
}

void ParMA_IO::errorPrint(const char* msg) {
    fprintf(stderr, "ERROR %s", msg);
}

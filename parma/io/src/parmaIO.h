#ifndef _PARMA_IO_H_
#define _PARMA_IO_H_

namespace ParMA_IO {
   void debugPrint(const char*, bool dbgOverride = false);
   void tracePrint(const char*, bool traceOverride = false);
   void statusPrint(const char*);
   void errorPrint(const char*);
};

#endif 

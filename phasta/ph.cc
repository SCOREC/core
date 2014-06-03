#include "ph.h"
#include <cstdio>
#include <cstdlib>
#include <stdarg.h>

namespace ph {

void fail(const char* format, ...)
{
  va_list ap;
  va_start(ap, format);
  vfprintf(stderr, format, ap);
  va_end(ap);
  fprintf(stderr,"\n");
  abort();
}

}

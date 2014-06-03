#include "ph.h"
#include <cstdio>
#include <cstdlib>
#include <stdarg.h>
#include <cassert>

/* both of these are for
   goToStepDir */
#include <sys/stat.h>
#include <unistd.h>

#include <sstream>

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

void goToStepDir(int step)
{
  std::stringstream ss;
  ss << step;
  std::string s = ss.str();
  int err = mkdir(s.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  assert(!err);
  err = chdir(s.c_str());
  assert(!err);
}

}

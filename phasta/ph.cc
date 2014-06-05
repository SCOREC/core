#include "ph.h"
#include <cstdio>
#include <cstdlib>
#include <stdarg.h>
#include <cassert>
#include <sstream>
#include <iostream>

/* OS-specific things try to stay here */
#include <sys/stat.h>
#include <unistd.h>

#include <PCU.h>

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

enum {
  DIR_MODE = S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH
};

void goToStepDir(int step)
{
  std::stringstream ss;
  ss << step;
  std::string s = ss.str();
  int err = mkdir(s.c_str(), DIR_MODE);
  assert(!err);
  err = chdir(s.c_str());
  assert(!err);
}

enum {
  DIR_FANOUT = 2048
};

std::string setupOutputDir(int parts)
{
  std::stringstream ss;
  ss << parts << "-procs_case";
  std::string s = ss.str();
  if (!PCU_Comm_Self())
    mkdir(s.c_str(), DIR_MODE);
  PCU_Barrier();
  return s;
}

void setupOutputSubdir(int parts, std::string& path)
{
  if (PCU_Comm_Peers() <= DIR_FANOUT)
    return;
  int self = PCU_Comm_Self();
  int subSelf = self % DIR_FANOUT;
  std::stringstream ss(path);
  ss << '/' << subSelf << '/';
  path = ss.str();
  if (!subSelf) {
    if (mkdir(path.c_str(), DIR_MODE)) {
      std::cerr << "overwriting directory " << path
          << " -Rank " << PCU_Comm_Self() << '\n';
    }
  }
  PCU_Barrier();
}

}

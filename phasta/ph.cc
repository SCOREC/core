#include "ph.h"
#include <cstdio>
#include <cstdlib>
#include <stdarg.h>
#include <cassert>
#include <sstream>
#include <iostream>
#include <fstream>

/* OS-specific things try to stay here */
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

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

static void my_mkdir(const char* name)
{
  int err = mkdir(name, DIR_MODE);
  if ((err == -1) && (errno == EEXIST)) {
    errno = 0;
    err = 0;
  }
  assert(!err);
}

static void my_chdir(const char* name)
{
  int err = chdir(name);
  assert(!err);
}

void goToStepDir(int step)
{
  std::stringstream ss;
  ss << step;
  std::string s = ss.str();
  if (!PCU_Comm_Self())
    my_mkdir(s.c_str());
  PCU_Barrier();
  my_chdir(s.c_str());
}

enum {
  DIR_FANOUT = 2048
};

std::string setupOutputDir()
{
  std::stringstream ss;
  ss << PCU_Comm_Peers() << "-procs_case/";
  std::string s = ss.str();
  if (!PCU_Comm_Self())
    my_mkdir(s.c_str());
  PCU_Barrier();
  return s;
}

void setupOutputSubdir(std::string& path)
{
  if (PCU_Comm_Peers() <= DIR_FANOUT)
    return;
  int self = PCU_Comm_Self();
  int subSelf = self % DIR_FANOUT;
  std::stringstream ss(path);
  ss << subSelf << '/';
  path = ss.str();
  if (!subSelf) {
    if (mkdir(path.c_str(), DIR_MODE)) {
      std::cerr << "overwriting directory " << path
          << " -Rank " << PCU_Comm_Self() << '\n';
    }
  }
  PCU_Barrier();
}

void writeAuxiliaryFiles(std::string path, int timestep)
{
  std::string numpePath = path;
  numpePath += "numpe.in";
  std::ofstream numpe(numpePath.c_str());
  assert(numpe.is_open());
  numpe << PCU_Comm_Peers() << '\n';
  numpe.close();
  std::string numstartPath = path;
  numstartPath += "numstart.dat";
  std::ofstream numstart(numstartPath.c_str());
  assert(numstart.is_open());
  numstart << timestep << '\n';
  numstart.close();
}

}

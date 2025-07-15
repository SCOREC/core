#include <lionPrint.h>

#include "ph.h"
#include <apfMDS.h>
#include <cstdio>
#include <cstdlib>
#include <stdarg.h>
#include <pcu_util.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstring>

#ifdef HAVE_SIMMETRIX
#include <phAttrib.h>
#include <apfSIM.h>
#include <gmi_sim.h>
#include <SimPartitionedMesh.h>
#endif

/* OS-specific things try to stay here */
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

namespace ph {

void fail(const char* format, ...)
{
  va_list ap;
  va_start(ap, format);
  lion_veprint(1, format, ap);
  va_end(ap);
  lion_eprint(1,"\n");
  abort();
}

enum {
  DIR_MODE = S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH
};

static bool my_mkdir(const char* name)
{
  errno = 0;
  int err = mkdir(name, DIR_MODE);
  if ((err == -1) && (errno == EEXIST)) {
    errno = 0;
    err = 0;
    return false;
  }
  PCU_ALWAYS_ASSERT(!err);
  return true;
}

static void my_chdir(const char* name)
{
  int err = chdir(name);
  PCU_ALWAYS_ASSERT(!err);
}

void goToParentDir() {
  my_chdir("..");
}

void goToStepDir(int step, pcu::PCU *PCUObj, bool all_mkdir)
{
  std::stringstream ss;
  ss << step;
  std::string s = ss.str();
  if (all_mkdir || !PCUObj->Self())
    my_mkdir(s.c_str());
  PCUObj->Barrier();
  my_chdir(s.c_str());
}

enum {
  DIR_FANOUT = 2048
};

std::string setupOutputDir(pcu::PCU *PCUObj, bool all_mkdir)
{
  std::stringstream ss;
  ss << PCUObj->Peers() << "-procs_case/";
  std::string s = ss.str();
  if (all_mkdir || !PCUObj->Self())
    my_mkdir(s.c_str());
  PCUObj->Barrier();
  return s;
}

void setupOutputSubdir(std::string& path, pcu::PCU *PCUObj, bool all_mkdir)
{
  if (PCUObj->Peers() <= DIR_FANOUT)
    return;
  int self = PCUObj->Self();
  int subSelf = self % DIR_FANOUT;
  int subGroup = self / DIR_FANOUT;
  std::stringstream ss;
  ss << path << subGroup << '/';
  path = ss.str();
  if (all_mkdir || !subSelf)
    my_mkdir(path.c_str());
  PCUObj->Barrier();
}

void setupInputSubdir(std::string& path, pcu::PCU *PCUObj)
{
  if (PCUObj->Peers() <= DIR_FANOUT)
    return;
  int self = PCUObj->Self();
  int subGroup = self / DIR_FANOUT;
  std::string newpath;
  std::stringstream ss;

  std::size_t found = path.find_last_of("/");
  if(found == std::string::npos) { 
    //no "/" character was found in path
    //this means no directory specified in path, only filename such as restart 
    //should not really happen since phasta files usually located in #-procs_case directories
    ss << "./" << subGroup << "/" << path;
  }
  else {
    //insert before the last "/" character the subgroup id 0, 1, 2 etc
    ss << path.substr(0,found) << "/" << subGroup << "/" << path.substr(found+1);
  }

//  lion_oprint(1,"Rank: %d - Path in setupInputSubdir: %s\n", self, ss.str().c_str());
  path = ss.str();
  PCUObj->Barrier();
}

void writeAuxiliaryFiles(std::string path, int timestep_or_dat, pcu::PCU *PCUObj)
{
  std::string numpePath = path;
  numpePath += "numpe.in";
  std::ofstream numpe(numpePath.c_str());
  PCU_ALWAYS_ASSERT(numpe.is_open());
  numpe << PCUObj->Peers() << '\n';
  numpe.close();
  std::string numstartPath = path;
  numstartPath += "numstart.dat";
  std::ofstream numstart(numstartPath.c_str());
  PCU_ALWAYS_ASSERT(numstart.is_open());
  numstart << timestep_or_dat << '\n';
  numstart.close();
}

bool mesh_has_ext(const char* filename, const char* ext)
{
  const char* c = strrchr(filename, '.');
  if (c) {
    ++c; /* exclude the dot itself */
    return !strcmp(c, ext);
  } else {
    return false;
  }
}

apf::Mesh2* loadMesh(gmi_model*& g, const char* meshfile, pcu::PCU *PCUObj) {
  apf::Mesh2* mesh;
#ifdef HAVE_SIMMETRIX
  /* if it is a simmetrix mesh */
  if (mesh_has_ext(meshfile, "sms")) {
    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    pGModel simModel = gmi_export_sim(g);
    pParMesh sim_mesh = PM_load(meshfile, simModel, progress);
    mesh = apf::createMesh(sim_mesh, PCUObj);

    Progress_delete(progress);
  } else
#endif
  /* if it is a SCOREC mesh */
  {
    mesh = apf::loadMdsMesh(g, meshfile, PCUObj);
  }
  return mesh;
}

}

#ifndef PH_H
#define PH_H

#include <apf.h>
#include <apfMesh2.h>
#include <string>

namespace ph {

void fail(const char* format, ...) __attribute__((noreturn,format(printf,1,2)));

void goToStepDir(int step, bool all_mkdir=false);
void goToParentDir();
std::string setupOutputDir(bool all_mkdir=false);
void setupInputSubdir(std::string& path);
void setupOutputSubdir(std::string& path, bool all_mkdir=false);
void writeAuxiliaryFiles(std::string path, int timestep_or_dat);
bool mesh_has_ext(const char* filename, const char* ext);
apf::Mesh2* loadMesh(gmi_model*& g, const char* meshfile);

}

#endif

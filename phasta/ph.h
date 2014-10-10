#ifndef PH_H
#define PH_H

#include <string>

namespace ph {

void fail(const char* format, ...) __attribute__((noreturn,format(printf,1,2)));

void goToStepDir(int step);
std::string setupOutputDir();
void setupInputSubdir(std::string& path);
void setupOutputSubdir(std::string& path);
void writeAuxiliaryFiles(std::string path, int timestep);

}

#endif

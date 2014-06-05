#ifndef PH_H
#define PH_H

#include <string>

namespace ph {

void fail(const char* format, ...) __attribute__((noreturn,format(printf,1,2)));

void goToStepDir(int step);
std::string setupOutputDir(int parts);
void setupOutputSubdir(int parts, std::string& path);

}

#endif

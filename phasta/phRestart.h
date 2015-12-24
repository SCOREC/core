#ifndef PH_RESTART_H
#define PH_RESTART_H

#include "phInput.h"
#include "phOutput.h"
#include <apfMesh.h>

namespace ph {

void readAndAttachSolution(Input& in, apf::Mesh* m);
void readAndAttachFields(Input& in, apf::Mesh* m);
void buildMapping(apf::Mesh* m);
void detachAndWriteSolution(Input& in, Output& out, 
    apf::Mesh* m, std::string path);
void attachZeroSolution(Input& in, apf::Mesh* m);

void detachField(apf::Field* f, double*& data, int& size);

}

#endif

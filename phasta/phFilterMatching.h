#ifndef PH_FILTER_MATCHING_H
#define PH_FILTER_MATCHING_H

#include "phBC.h"
#include "phInput.h"
#include <apfMesh2.h>

namespace ph {

void enterFilteredMatching(apf::Mesh2* m, Input& in, BCs& bcs);
void exitFilteredMatching(apf::Mesh2* m);

}

#endif

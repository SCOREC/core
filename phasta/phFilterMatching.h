#ifndef PH_FILTER_MATCHING_H
#define PH_FILTER_MATCHING_H

#include "phBC.h"
#include <apfMesh2.h>
#include <vector>

namespace ph {

typedef std::vector<apf::Matches> SavedMatches;

void saveMatches(apf::Mesh* m, int dim, SavedMatches& sm);
void restoreMatches(apf::Mesh2* m, int dim, SavedMatches& sm);
void filterMatches(apf::Mesh2* m, BCs& bcs);

}

#endif

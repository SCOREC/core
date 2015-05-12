#ifndef PH_FILTER_MATCHING_H
#define PH_FILTER_MATCHING_H

#include "phBC.h"
#include <apfMesh2.h>
#include <vector>

namespace ph {

typedef std::vector<apf::Matches> SavedMatches;
typedef std::set<gmi_ent*> ModelSet;
typedef std::map<gmi_ent*, ModelSet> ModelMatching;

void saveMatches(apf::Mesh* m, int dim, SavedMatches& sm);
void restoreMatches(apf::Mesh2* m, int dim, SavedMatches& sm);
void getFullAttributeMatching(gmi_model* m, BCs& bcs, ModelMatching& mm);
void filterMatching(apf::Mesh2* m, ModelMatching& mm, int dim);

}

#endif

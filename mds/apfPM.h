/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#ifndef MDS_PART_H
#define MDS_PART_H

#include <apfMesh2.h>

namespace apf {

struct PME
{
  PME(Parts const& i)
  {
    owner = -1;
    ids.assign(i.begin(), i.end());
    refs = 0;
  }
  bool operator<(PME const& other) const
  {
    return ids < other.ids;
  }
  int owner;
  std::vector<int> ids;
  int refs;
};

typedef std::set<PME> PM;

PME* getPME(PM& ps, Parts const& ids);
void putPME(PM& ps, PME* p);
void updateOwners(Mesh* m, PM& ps);

void initResidence(Mesh2* m, int dim);
void stitchMesh(Mesh2* m);

}

#endif

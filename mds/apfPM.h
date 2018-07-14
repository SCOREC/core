/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#ifndef APF_PM_H
#define APF_PM_H

#include <apfMesh2.h>

namespace apf {

struct PME
{
  PME(int n, Parts const& i, int o)
  {
    ID=n;
    owner=-1;
    if (o>-1) owner = o;
    ids.assign(i.begin(), i.end());
    refs = 0;
#ifdef DEBUG
  printf("(%d) create a new PME %d (pids.size=%d, ownder=%d)\n", PCU_Comm_Self(), ID, ids.size(), owner);
#endif
  }

  bool operator<(PME const& other) const
  {
    return ids < other.ids;
  }
  int owner;
  std::vector<int> ids;
  int refs;
  int ID;
};

typedef std::set<PME> PM;

void deletePM (PM& ps);
void deletePMent(PM& ps, PME* p);
PME* getPMent(PM& ps, apf::Parts const& pids, int);

PME* getPME(PM& ps, Parts const& ids);
void putPME(PM& ps, PME* p);
void updateOwners(Mesh* m, PM& ps);

void stitchMesh(Mesh2* m);

void remapPM(PM& pm,
    int (*map)(int, void*), void* user);

}

#endif

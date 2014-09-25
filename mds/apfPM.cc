/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#include <PCU.h>
#include "apfPM.h"
#include <apf.h>

namespace apf {

typedef std::map<int,size_t> CountMap;

PME* getPME(PM& ps, apf::Parts const& ids)
{
  PME const& cp = *(ps.insert(PME(ids)).first);
  /* always annoyed by this flaw in std::set */
  PME& p = const_cast<PME&>(cp);
  ++(p.refs);
  return &p;
}

void putPME(PM& ps, PME* p)
{
  --(p->refs);
  if (!(p->refs))
    ps.erase(*p);
}

static void getAdjacentParts(apf::Mesh* m, PM& ps, apf::Parts& ids)
{
  APF_ITERATE(PM, ps, it) {
    std::vector<int> const& rp = it->ids;
    ids.insert(rp.begin(), rp.end());
  }
  ids.erase(m->getId());
}

static void getCountMap(apf::Mesh* m, PM& ps, CountMap& mp)
{
  apf::Parts peers;
  size_t n;
  getAdjacentParts(m, ps, peers);
  n = m->count(m->getDimension());
  PCU_Comm_Begin();
  APF_ITERATE(apf::Parts, peers, it)
    PCU_COMM_PACK(*it, n);
  PCU_Comm_Send();
  mp[m->getId()] = n;
  while (PCU_Comm_Listen()) {
    PCU_COMM_UNPACK(n);
    mp[PCU_Comm_Sender()] = n;
  }
}

static void setOwners(PM& ps, CountMap& mp)
{
  APF_ITERATE(PM, ps, it) {
    PME const& cp = *it;
    PME& p = const_cast<PME&>(cp); /* again with the silly */
    std::vector<int> const& ids = p.ids;
    assert(ids.size());
    int owner = ids[0];
    assert(mp.count(owner));
    for (size_t i = 1; i < ids.size(); ++i) {
      assert(mp.count(ids[i]));
      if (mp[ids[i]] < mp[owner])
        owner = ids[i];
    }
    p.owner = owner;
  }
}

void updateOwners(apf::Mesh* m, PM& ps)
{
  CountMap mp;
  getCountMap(m, ps, mp);
  setOwners(ps, mp);
}

void remapPM(PM& pm,
    int (*map)(int, void*), void* user)
{
  APF_ITERATE(PM, pm, it) {
    PME const& cp = *it;
    PME& p = const_cast<PME&>(cp); /* yep */
    p.owner = map(p.owner, user);
    std::vector<int>& ids = p.ids;
/* note: we can only do this because operator<(std::vector<T>...)
   uses lexicographical comparison, and so for vectors A and B,
   A < B does not change if all the elements of A and B are
   multiplied or divided by a constant factor, so long as
   the resulting ids are also unique.
   
   any map which changes the results of lexicographical
   comparison breaks the ordering of PME's in the PM. */
    for (size_t i = 0; i < ids.size(); ++i)
      ids[i] = map(ids[i], user);
  }
}

}

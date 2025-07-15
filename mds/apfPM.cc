/****************************************************************************** 

  Copyright 2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/


#include "apfPM.h"
#include <apf.h>
#include <pcu_util.h>

namespace apf {

// the following three functions are for exclusive use by pumi_mesh_loadAll
// the last arg "owner" is used only if a new pmodel entity is created
PME* getPMent(PM& ps, apf::Parts const& pids, int owner)
{
  APF_ITERATE(PM, ps, it) 
  {
    bool equal=true;
    APF_ITERATE(Parts,pids,pit)
    {
      bool found=false;
      for (size_t i = 0; i < (*it).ids.size(); ++i)
        if ((*it).ids[i]==*pit)
          found=true;
      if (!found)
      {
        equal=false;
        break;
      }
    }

    if (equal)
    {  
      PME& p = const_cast<PME&>(*it);
      ++(p.refs);
      return &p;
    }
  }
  static int pme_id=ps.size();
  PME *pme = new PME(pme_id++, pids, owner);
  ps.insert(*pme);
  ++(pme->refs);
  return pme;
}

// the partition classification of mesh entities has to be updated separately
void deletePM(PM& ps)
{  
  APF_ITERATE(PM, ps, it) 
    ps.erase(*it);
}

void deletePMent(PM& ps, PME* p)
{
  ps.erase(*p);
}


typedef std::map<int,size_t> CountMap;

PME* getPME(PM& ps, apf::Parts const& ids)
{
  PME const& cp = *(ps.insert(PME(ps.size(), ids, -1)).first);
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
  m->getPCU()->Begin();
  APF_ITERATE(apf::Parts, peers, it)
    m->getPCU()->Pack(*it, n);
  m->getPCU()->Send();
  mp[m->getId()] = n;
  while (m->getPCU()->Listen()) {
    m->getPCU()->Unpack(n);
    mp[m->getPCU()->Sender()] = n;
  }
}

static void setOwners(PM& ps, CountMap& mp)
{
  APF_ITERATE(PM, ps, it) {
    PME const& cp = *it;
    PME& p = const_cast<PME&>(cp); /* again with the silly */
    std::vector<int> const& ids = p.ids;
    PCU_ALWAYS_ASSERT(ids.size());
    int owner = ids[0];
//    PCU_ALWAYS_ASSERT(mp.count(owner)); // seol - this doesn't work for ghost copy
    for (size_t i = 1; i < ids.size(); ++i) 
    {
//      PCU_ALWAYS_ASSERT(mp.count(ids[i]));  // seol - this doesn't work for ghost copy
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

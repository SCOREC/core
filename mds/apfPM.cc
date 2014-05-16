/****************************************************************************** 

  (c) 2011-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#include "apfPM.h"
#include <apf.h>
#include <PCU.h>

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
    for (size_t i = 0; i < ids.size(); ++i) {
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

void initResidence(apf::Mesh2* m, int d)
{
  apf::MeshIterator* it = m->begin(d);
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    apf::Copies remotes;
    m->getRemotes(e, remotes);
    apf::Parts parts;
    APF_ITERATE(apf::Copies, remotes, rit)
      parts.insert(rit->first);
    parts.insert(m->getId());
    m->setResidence(e, parts);
  }
  m->end(it);
}

static void intersect(
    Parts& a,
    Parts const& b)
{
  for (Parts::iterator it = a.begin(); it != a.end();) {
    if ( ! b.count(*it))
      a.erase(*(it++));
    else
      ++it;
  }
}

static void getCandidateParts(Mesh* m, MeshEntity* e, Parts& parts)
{
  int d = getDimension(m, e);
  Downward down;
  int nd = m->getDownward(e, d - 1, down);
  //down adjacencies have to at least be shared
  for (int i = 0; i < nd; ++i)
    if (!m->isShared(down[i]))
      return;
  bool first = true;
  for (int i = 0; i < nd; ++i) {
    MeshEntity* da = down[i];
    Parts rp;
    m->getResidence(da, rp);
    if (first) {
      parts = rp;
      first = false;
    } else {
      intersect(parts,rp);
    }
  }
}

static void packProposal(Mesh* m, MeshEntity* e, int to)
{
  int t = m->getType(e);
  PCU_COMM_PACK(to,t);
  PCU_COMM_PACK(to,e);
  int d = getDimension(m, e);
  Downward down;
  int nd = m->getDownward(e, d - 1, down);
  PCU_COMM_PACK(to,nd);
  for (int i = 0; i < nd; ++i) {
    Copies remotes;
    m->getRemotes(down[i], remotes);
    MeshEntity* dr = remotes[to];
    PCU_COMM_PACK(to,dr);
  }
}

static void unpackProposal(int& t, MeshEntity*& e, Downward& da)
{
  PCU_COMM_UNPACK(t);
  PCU_COMM_UNPACK(e);
  int nd;
  PCU_COMM_UNPACK(nd);
  for (int i=0; i < nd; ++i)
    PCU_COMM_UNPACK(da[i]);
}

void stitchMesh(Mesh2* m)
{
  initResidence(m, 0);
  int d_max = m->getDimension();
  MeshEntity* e;
  for (int d=1; d < d_max; ++d) {
    PCU_Comm_Begin();
    MeshIterator* it = m->begin(d);
    while ((e = m->iterate(it))) {
      Parts candidateParts;
      getCandidateParts(m, e, candidateParts);
      candidateParts.erase(m->getId());
      APF_ITERATE(Parts, candidateParts, pit)
        packProposal(m, e,*pit);
    }
    m->end(it);
    PCU_Comm_Send();
    while (PCU_Comm_Listen()) {
      int from = PCU_Comm_Sender();
      while (!PCU_Comm_Unpacked()) {
        int t;
        Downward da;
        unpackProposal(t, e, da);
        MeshEntity* found = findUpward(m, t, da);
        if (found)
          m->addRemote(found, from, e);
      }
    }
    initResidence(m, d);
  }
}

}

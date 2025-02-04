/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maMatch.h"
#include "maMesh.h"
#include "maAdapt.h"
#include "maRefine.h"

namespace ma {

static void packSplits(int to, EntityArray& splits, pcu::PCU *PCUObj)
{
  PCUObj->Pack(to,
      static_cast<void*>(&(splits[0])),
      splits.getSize()*sizeof(Entity*));
}

static void unpackSplits(EntityArray& splits, pcu::PCU *PCUObj)
{
  PCUObj->Unpack(
      static_cast<void*>(&(splits[0])),
      splits.getSize()*sizeof(Entity*));
}

void matchNewElements(Refine* r)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  long face_count = 0;
  for (int d=1; d < m->getDimension(); ++d)
  {
    m->getPCU()->Begin();
    for (size_t i=0; i < r->toSplit[d].getSize(); ++i)
    {
      Entity* e = r->toSplit[d][i];
      apf::Matches matches;
      m->getMatches(e,matches);
      if ( ! matches.getSize())
        continue;
      EntityArray& splits = r->newEntities[d][i];
      for (size_t i=0; i < matches.getSize(); ++i)
      {
        int to = matches[i].peer;
        Entity* match = matches[i].entity;
        m->getPCU()->Pack(to,match);
        packSplits(to,splits,m->getPCU());
      }
    }
    m->getPCU()->Send();
    while (m->getPCU()->Listen())
    {
      int from = m->getPCU()->Sender();
      while ( ! m->getPCU()->Unpacked())
      {
        Entity* e;
        m->getPCU()->Unpack(e);
        int number;
        m->getIntTag(e,r->numberTag,&number);
        EntityArray& splits = r->newEntities[d][number];
        EntityArray remoteSplits(splits.getSize());
        unpackSplits(remoteSplits,m->getPCU());
        for (size_t i=0; i < splits.getSize(); ++i)
          m->addMatch(splits[i],from,remoteSplits[i]);
        if (d==2) ++face_count;
      }
    }
  }
  face_count = m->getPCU()->Add<long>(face_count);
  print(m->getPCU(), "updated matching for %li faces", face_count);
}

void preventMatchedCavityMods(Adapt* a)
{
  Mesh* m = a->mesh;
  if (!m->hasMatching())
    return;
  Iterator* it = m->begin(0);
  Entity* v;
  while ((v = m->iterate(it))) {
    apf::Matches matches;
    m->getMatches(v, matches);
    if (!matches.getSize())
      continue;
    apf::Up edges;
    setFlag(a, v, DONT_COLLAPSE);
    m->getUp(v, edges);
    for (int i = 0; i < edges.n; ++i)
      setFlag(a, edges.e[i], DONT_COLLAPSE | DONT_SWAP);
  }
  m->end(it);
}

}

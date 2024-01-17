/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maMatchedSnapper.h"
#include "maAdapt.h"
#include "maShapeHandler.h"
#include <apfCavityOp.h>
#include <pcu_util.h>
#include <PCU.h>
#include <iostream>

namespace ma {

MatchedSnapper::MatchedSnapper(Adapt* a, Tag* st, bool is)
{
  adapter = a;
  snapTag = st;
  isSimple = is;
  sharing = apf::getSharing(a->mesh);
  vert = 0;
}

MatchedSnapper::~MatchedSnapper()
{
  delete sharing;
  // delete snapper objects in snappers
  for (size_t i = 0; i < snappers.getSize(); i++) {
    delete snappers[i];
  }
}

void MatchedSnapper::setVert(Entity* v)
{
  snappers.setSize(0);
  snappers.setSize(1);
  snappers[0] = new Snapper(adapter, snapTag, isSimple);
  snappers[0]->setVert(v);
  vert = v;
}

bool MatchedSnapper::requestLocality(apf::CavityOp* o)
{
  return snappers[0]->requestLocality(o);
}

void MatchedSnapper::setVerts()
{
  Entity* v = vert;
  delete snappers[0];
  snappers.setSize(0);

  apf::CopyArray copies;
  sharing->getCopies(v, copies);
  APF_ITERATE(apf::CopyArray, copies, it) {
    PCU_ALWAYS_ASSERT(it->peer == adapter->mesh->getPCU()->Self());
    PCU_ALWAYS_ASSERT(it->entity != v);
  }

  snappers.setSize(copies.getSize() + 1);
  for (unsigned i = 0; i < snappers.getSize(); i++)
    snappers[i] = new Snapper(adapter, snapTag, isSimple);
  snappers[0]->setVert(v);
  for (unsigned i = 0; i < copies.getSize(); i++)
    snappers[i+1]->setVert(copies[i].entity);

  // cache the original locations
  locations.setSize(snappers.getSize());
  locations[0] = getPosition(adapter->mesh, v);
  for (unsigned i = 0; i < copies.getSize(); i++)
    locations[i+1] = getPosition(adapter->mesh, copies[i].entity);
}

bool MatchedSnapper::trySnaps()
{
  // if only 1, use the old snappers run() method
  if (snappers.getSize() == 1)
    return snappers[0]->run();
  bool allSnapped = true;
  for (unsigned i = 0; i < snappers.getSize(); i++) {
    allSnapped = allSnapped && snappers[i]->trySimpleSnap();
    if (! allSnapped)
      break;
  }
  if (! allSnapped) {
    cancelSnaps();
    return false;
  }
  return true;
}

void MatchedSnapper::cancelSnaps()
{
  Mesh* m = adapter->mesh;
  for (unsigned i = 0; i < snappers.getSize(); i++)
    m->setPoint(snappers[i]->getVert(), 0, locations[i]);
}

}

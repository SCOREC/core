/****************************************************************************** 

  Copyright 2015 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maMatchedCollapse.h"
#include "maAdapt.h"
#include <cassert>
#include <PCU.h>

namespace ma {

MatchedCollapse::MatchedCollapse(Adapt* a):
  adapt(a)
{
  mesh = a->mesh;
  sharing = apf::getSharing(mesh);
}

void MatchedCollapse::setEdge(Entity* e)
{
  collapses.setSize(0);
  collapses.setSize(1);
  collapses[0].Init(adapt);
  bool ok = collapses[0].setEdge(e);
  assert(ok);
}

bool MatchedCollapse::requestLocality(apf::CavityOp* o)
{
  return collapses[0].requestLocality(o);
}

void MatchedCollapse::setEdges()
{
  Entity* e = collapses[0].edge;
  collapses.setSize(0);
  apf::CopyArray copies;
  sharing->getCopies(e, copies);
  APF_ITERATE(apf::CopyArray, copies, it) {
    assert(it->peer == PCU_Comm_Self());
    assert(it->entity != e);
  }
  collapses.setSize(copies.getSize() + 1);
  for (unsigned i = 0; i < collapses.getSize(); ++i)
    collapses[i].Init(adapt);
  bool ok = collapses[0].setEdge(e);
  assert(ok);
  for (unsigned i = 0; i < copies.getSize(); ++i) {
    ok = collapses[i + 1].setEdge(copies[i].entity);
    assert(ok);
  }
}

bool MatchedCollapse::checkTopo()
{
  bool ok = true;
  for (unsigned i = 0; i < collapses.getSize(); ++i)
    if (!collapses[i].checkTopo())
      ok = false;
  if (!ok)
    for (unsigned i = 0; i < collapses.getSize(); ++i)
      collapses[i].unmark();
  return ok;
}

/*
bool MatchedCollapse::tryBothDirections(double qualityToBeat)
{
  computeElementSets();
  if (tryThisDirection(qualityToBeat))
    return true;
  for (unsigned i = 0; i < collapses.getSize(); ++i)
    if (!getFlag(adapt, collapses[i].vertToKeep, COLLAPSE))
      return false;
  for (unsigned i = 0; i < collapses.getSize(); ++i)
    std::swap(collapses[i].vertToKeep, collapses[i].vertToCollapse);
  computeElementSets();
  return tryThisDirection(qualityToBeat);
}
*/

}

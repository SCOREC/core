/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maDoubleSplitCollapse.h"
#include "maAdapt.h"
#include "maShape.h"

namespace ma {

DoubleSplitCollapse::DoubleSplitCollapse(Adapt* a):
  splits(a)
{
  collapse.Init(a);
}

void DoubleSplitCollapse::getNewElements(EntityArray& e)
{
  Mesh* m = getAdapt()->mesh;
  EntityArray& c = collapse.newElements;
  EntitySet& b = collapse.elementsToCollapse;
  EntityArray ab;
  m->getAdjacent(collapse.vertToKeep,3,ab);
  e.setSize(ab.getSize() - b.size() + c.getSize());
  size_t i=0;
  for (size_t j=0; j < ab.getSize(); ++j)
    if ( ! b.count(ab[j]))
      e[i++] = ab[j];
  for (size_t j=0; j < c.getSize(); ++j)
    e[i++] = c[j];
  assert(i==e.getSize());
}

bool DoubleSplitCollapse::tryBothCollapses(Entity* e)
{
  if ( ! collapse.setEdge(e))
    return false;
  if ( ! collapse.checkClass())
    return false;
  if ( ! collapse.checkTopo())
    return false;
  return collapse.tryBothDirections(oldQuality);
}

void DoubleSplitCollapse::accept()
{
  splits.destroyOldElements();
  collapse.destroyOldElements();
}

bool DoubleSplitCollapse::run(Entity** edges)
{
  Adapt* a = getAdapt();
  Mesh* m = a->mesh;
  if ( ! splits.setEdges(edges,2))
    return false;
  oldQuality = getWorstQuality(a,splits.getTets());
  splits.makeNewElements();
  splits.transfer();
  Entity* splitVerts[2];
  splitVerts[0] = splits.getSplitVert(0);
  splitVerts[1] = splits.getSplitVert(1);
  Entity* edge = findUpward(m, apf::Mesh::EDGE, splitVerts);
  if (tryBothCollapses(edge))
  {
    accept();
    return true;
  }
  else
  {
    splits.cancel();
    return false;
  }
}

Adapt* DoubleSplitCollapse::getAdapt()
{
  return collapse.adapt;
}

}

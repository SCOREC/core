/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maDoubleSplitCollapse.h"
#include "maAdapt.h"
#include "maShape.h"
#include <pcu_util.h>

namespace ma {

DoubleSplitCollapse::DoubleSplitCollapse(Adapt* a):
  splits(a)
{
  collapse.Init(a);
  oldQuality = 2;
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
  PCU_ALWAYS_ASSERT(i==e.getSize());
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

  /* This is to make sure the quality of newly created elements
   * by the double splits is not less than the in->goodQuality
   */
  if (a->input->shouldCheckQualityForDoubleSplits) {
    EntityArray toCheck;
    Upward regionsToSplitVert0;
    Upward regionsToSplitVert1;
    Upward regionsToEdge;
    m->getAdjacent(splitVerts[0], m->getDimension(), regionsToSplitVert0);
    m->getAdjacent(splitVerts[1], m->getDimension(), regionsToSplitVert1);
    m->getAdjacent(edge, m->getDimension(), regionsToEdge);
    PCU_ALWAYS_ASSERT(regionsToEdge.getSize() == 4);
    toCheck.setSize(regionsToSplitVert0.getSize() +
		    regionsToSplitVert1.getSize() -
		    2*regionsToEdge.getSize());
    size_t count = 0;
    for (size_t i = 0; i < regionsToSplitVert0.getSize(); i++) {
      Entity* region = regionsToSplitVert0[i];
      if (region == regionsToEdge[0] ||
	  region == regionsToEdge[1] ||
	  region == regionsToEdge[2] ||
	  region == regionsToEdge[3])
	continue;
      else {
	toCheck[count] = region;
	count++;
      }
    }
    for (size_t i = 0; i < regionsToSplitVert1.getSize(); i++) {
      Entity* region = regionsToSplitVert1[i];
      if (region == regionsToEdge[0] ||
	  region == regionsToEdge[1] ||
	  region == regionsToEdge[2] ||
	  region == regionsToEdge[3])
	continue;
      else {
	toCheck[count] = region;
	count++;
      }
    }
    PCU_ALWAYS_ASSERT(count == toCheck.getSize());
    double toCheckQuality = getWorstQuality(a, toCheck);
    if (toCheckQuality < a->input->goodQuality) {
      splits.cancel();
      return false;
    }
  }

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

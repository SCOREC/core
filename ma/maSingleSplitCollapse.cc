/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maSingleSplitCollapse.h"
#include "maAdapt.h"
#include "maShape.h"
#include <pcu_util.h>

namespace ma {

SingleSplitCollapse::SingleSplitCollapse(Adapt* a):
  splits(a)
{
  collapse.Init(a);
  oldQuality = 2;
}

void SingleSplitCollapse::getNewElements(EntityArray& e)
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

bool SingleSplitCollapse::tryBothCollapses(Entity* e)
{
  if ( ! collapse.setEdge(e))
    return false;
  if ( ! collapse.checkClass())
    return false;
  if ( ! collapse.checkTopo())
    return false;
  EntityArray& preSplit = splits.getTets();
  for (size_t i=0; i < preSplit.getSize(); ++i)
    collapse.elementsToIgnore.insert(preSplit[i]);
  return collapse.tryBothDirections(oldQuality);
  collapse.elementsToIgnore.clear();
}

void SingleSplitCollapse::accept()
{
  splits.destroyOldElements();
  collapse.destroyOldElements();
}

bool SingleSplitCollapse::run(Entity* edge, Entity* vert)
{
  oldEdge = edge;
  oldVert = vert;
  Adapt* a = getAdapt();
  Mesh* m = a->mesh;
  if ( ! splits.setEdges(&edge,1))
    return false;
  oldQuality = getWorstQuality(a,splits.getTets());
  splits.makeNewElements();
  splits.transfer();
  Entity* collVerts[2];
  collVerts[0] = splits.getSplitVert(0);
  collVerts[1] = vert;
  Entity* collEdge = findUpward(m, apf::Mesh::EDGE, collVerts);

  /* This is to make sure the quality of newly created elements
   * by the edge split is not less than in->goodQuality
   */
  if (a->input->shouldCheckQualityForDoubleSplits) {
    EntityArray toCheck;
    Upward regionsToCollVert0;
    Upward regionsToEdge;
    m->getAdjacent(collVerts[0], m->getDimension(), regionsToCollVert0);
    m->getAdjacent(collEdge, m->getDimension(), regionsToEdge);
    // there would be 2 regions if on the boundary, 4 if interior
    PCU_ALWAYS_ASSERT(regionsToEdge.getSize() == 2 ||
                      regionsToEdge.getSize() == 4);
    toCheck.setSize(regionsToCollVert0.getSize() -
		    regionsToEdge.getSize());
    size_t count = 0;
    for (size_t i = 0; i < regionsToCollVert0.getSize(); i++) {
      Entity* region = regionsToCollVert0[i];
      // ignoring these because they will be collapsed
      if (region == regionsToEdge[0] ||
	  region == regionsToEdge[1] ||
          (regionsToEdge.getSize() == 4 &&
           (region == regionsToEdge[2] ||
            region == regionsToEdge[3])))
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

  if (tryBothCollapses(collEdge))
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

Adapt* SingleSplitCollapse::getAdapt()
{
  return collapse.adapt;
}

void SingleSplitCollapse::IgnoringCollapse::computeElementSets()
{
  Upward adjacent;
  Mesh* m = adapt->mesh;
  m->getAdjacent(edge,m->getDimension(),adjacent);
  elementsToCollapse.clear();
  APF_ITERATE(Upward,adjacent,it)
    if ( ! elementsToIgnore.count(*it))
      elementsToCollapse.insert(*it);
  m->getAdjacent(vertToCollapse,m->getDimension(),adjacent);
  elementsToKeep.clear();
  APF_ITERATE(Upward,adjacent,it)
    if ( ! (elementsToCollapse.count(*it) ||
            elementsToIgnore.count(*it)))
      elementsToKeep.insert(*it);
  PCU_ALWAYS_ASSERT(elementsToKeep.size());
}

bool SingleSplitCollapse::IgnoringCollapse::setEdge(Entity* e)
{
  if (getFlag(adapt,e,DONT_COLLAPSE))
    return false;
  edge = e;
  vertToCollapse = 0;
  vertToKeep = 0;
  elementsToCollapse.clear();
  elementsToKeep.clear();
  elementsToIgnore.clear();
  return true;
}

}

/****************************************************************************** 

  Copyright (c) 2004-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maVertRemover.h"
#include "maAdapt.h"
#include "maShape.h"

namespace ma {

void VertRemover::Init(Adapt* a)
{
  adapter = a;
  mesh = a->mesh;
  collapse.Init(a);
}

void VertRemover::setVert(Entity* v)
{
  vert = v;
}

void VertRemover::findEdges()
{
  mesh->getAdjacent(vert,1,edges);
}

bool VertRemover::didImproveQuality()
{
  EntityArray oldElements;
  collapse.getOldElements(oldElements);
  EntityArray& newElements = collapse.newElements;
  double oldQuality = getWorstQuality(adapter,oldElements);
  double newQuality = getWorstQuality(adapter,newElements);
  return newQuality > oldQuality;
}

bool VertRemover::tryToCollapse(Entity* e)
{
  if ( ! collapse.setEdge(e))
    return false;
  if ( ! collapse.checkClass())
    return false;
  if ( ! collapse.checkTopo())
    return false;
  if ( ! getFlag(adapter,vert,COLLAPSE))
  {
    collapse.unmark();
    return false;
  }
  if (collapse.vertToCollapse != vert)
    std::swap(collapse.vertToCollapse,collapse.vertToKeep);
  assert(collapse.vertToCollapse==vert);
  double oldQuality = collapse.getOldQuality();
  if ( ! collapse.tryThisDirection(oldQuality))
    return false;
  collapse.destroyOldElements();
  return true;
}

bool VertRemover::run()
{
  for (size_t i=0; i < edges.getSize(); ++i)
    if (tryToCollapse(edges[i]))
      return true;
  return false;
}

}

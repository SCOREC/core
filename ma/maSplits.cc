/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maSplits.h"
#include "maAdapt.h"

namespace ma {

Splits::Splits(Adapt* a):
  refiner(a->refine)
{
  resetCollection(refiner);
  collectForTransfer(refiner);
}

bool Splits::setEdges(Entity** e, int edgeCount)
{
  Adapt* a = getAdapt();
  for (int i=0; i < edgeCount; ++i)
    if (getFlag(a,e[i],DONT_SPLIT))
      return false;
  for (int d=1; d <= 3; ++d)
    refiner->toSplit[d].setSize(0);
  int n[4] = {0,0,0,0};
  for (int i=0; i < edgeCount; ++i)
    addEdgePreAllocation(refiner,e[i],n);
  allocateRefine(refiner,n);
  n[0]=n[1]=n[2]=n[3]=0;
  for (int i=0; i < edgeCount; ++i)
  {
    setFlag(a,e[i],SPLIT);
    addEdgePostAllocation(refiner,e[i],n);
  }
  return true;
}

void Splits::makeNewElements()
{
  splitElements(refiner);
}

void Splits::cancel()
{
  Adapt* a = getAdapt();
  Mesh* m = a->mesh;
  EntitySet deleting;
  EntityArray& edges = refiner->toSplit[1];
  for (size_t i=0; i < edges.getSize(); ++i)
  {
    Upward elements;
    m->getAdjacent(getSplitVert(i),m->getDimension(),elements);
    for (size_t j=0; j < elements.getSize(); ++j)
      deleting.insert(elements[j]);
  }
  APF_ITERATE(EntitySet,deleting,it)
    destroyElement(a,*it);
  for (int d=1; d <= m->getDimension(); ++d)
  {
    EntityArray& e = refiner->toSplit[d];
    for (size_t i=0; i < e.getSize(); ++i)
    {
      clearFlag(a,e[i],SPLIT);
      m->removeTag(e[i],refiner->numberTag);
    }
  }
  forgetNewEntities(refiner);
}

void Splits::transfer()
{
  transferElements(refiner);
}

void Splits::destroyOldElements()
{
  destroySplitElements(refiner);
  cleanSplitVerts(refiner);
  forgetNewEntities(refiner);
}

Entity* Splits::getSplitVert(int i)
{
  return findSplitVert(refiner,1,i);
}

}

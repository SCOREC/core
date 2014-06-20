/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maVertRemover.h"
#include "maAdapt.h"

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
  mesh->getUp(vert,edges);
}

bool VertRemover::tryToCollapse(Entity* e)
{
  if (!setupCollapse(collapse, e, vert))
    return false;
  double oldQuality = collapse.getOldQuality();
  if ( ! collapse.tryThisDirection(oldQuality))
    return false;
  collapse.destroyOldElements();
  return true;
}

bool VertRemover::run()
{
  for (int i = 0; i < edges.n; ++i)
    if (tryToCollapse(edges.e[i]))
      return true;
  return false;
}

}

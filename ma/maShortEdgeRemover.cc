/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maShortEdgeRemover.h"
#include "maShape.h"
#include "maAdapt.h"
#include <apfCavityOp.h>

namespace ma {

ShortEdgeRemover::ShortEdgeRemover(Adapt* a)
{
  adapter = a;
  mesh = a->mesh;
  for (int i=0; i < 2; ++i)
    vertRemovers[i].Init(a);
}

void ShortEdgeRemover::setEdge(Entity* e)
{
  edge = e;
  Entity* verts[2];
  mesh->getDownward(edge,0,verts);
  vertRemovers[0].setVert(verts[0]);
  vertRemovers[1].setVert(verts[1]);
}

void ShortEdgeRemover::findEdges()
{
  vertRemovers[0].findEdges();
  vertRemovers[1].findEdges();
}

bool ShortEdgeRemover::requestLocality(apf::CavityOp* o)
{
  Entity* verts[2] = {vertRemovers[0].getVert(),
                      vertRemovers[1].getVert()};
  if (mesh->isShared(verts[0])||mesh->isShared(verts[1]))
    return o->requestLocality(verts,2);
  findEdges();
  apf::Up* edges[2] = { & (vertRemovers[0].getEdges()),
                        & (vertRemovers[1].getEdges()) };
  EntityArray otherVerts(edges[0]->n + edges[1]->n - 2);
  size_t k=0;
  for (int i=0; i < 2; ++i)
    for (int j=0; j < edges[i]->n; ++j)
    {
      Entity* e = edges[i]->e[j];
      if (e == edge) continue;
      otherVerts[k++] = getEdgeVertOppositeVert(mesh,e,verts[i]);
    }
  assert(k==otherVerts.getSize());
/* omitting the two original verts does not remove elements
   from the cavity being requested */
  return o->requestLocality(&(otherVerts[0]),otherVerts.getSize());
}

bool ShortEdgeRemover::tryToRemoveVert(int vi)
{
  return vertRemovers[vi].run();
}

bool ShortEdgeRemover::run()
{
  for (int i=0; i < 2; ++i)
    if (tryToRemoveVert(i))
      return true;
  return false;
}

}

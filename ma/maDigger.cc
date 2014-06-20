/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maDigger.h"
#include "maAdapt.h"
#include <apfCavityOp.h>
#include <map>

namespace ma {

Digger::Digger(Adapt* a, Tag* st)
{
  adapter = a;
  mesh = a->mesh;
  snapTag = st;
  collapse.Init(a);
}

bool Digger::setVert(Entity* v, apf::CavityOp* o)
{
  vert = v;
  if (!o->requestLocality(&vert, 1))
    return false;
  mesh->getUp(vert,edges);
  return o->requestLocality(&edges.e[0], edges.n);
}

bool Digger::tryToCollapse(Entity* e)
{
  Entity* ov = apf::getEdgeVertOppositeVert(mesh, e, vert);
  if (!setupCollapse(collapse, e, ov))
    return false;
  double q = adapter->input->validQuality;
  if ( ! collapse.tryThisDirection(q))
    return false;
  collapse.destroyOldElements();
  return true;
}

/* instead of trying edges at random, let's prefer edges
   which are going in the direction in which we are snapping.
   For that, we can take the dot product of the vector along
   the edge and the intended snap vector. */

typedef std::map<double, Entity*, std::greater<double> > Map;
/* put the most aligned edge first ----^ */

static void sortEdges(
    Mesh* mesh,
    Tag* snapTag,
    Entity* vert,
    apf::Up& edges,
    Map& map)
{
  Vector x = getPosition(mesh, vert);
  Vector s;
  mesh->getDoubleTag(vert, snapTag, &s[0]);
  Vector snapVector = s - x;
  for (int i = 0; i < edges.n; ++i) {
    Entity* e = edges.e[i];
    Entity* ov = apf::getEdgeVertOppositeVert(mesh, e, vert);
    Vector ox = getPosition(mesh, ov);
    Vector edgeVector = ox - x;
    double d = edgeVector * snapVector;
    if (d > 0) /* ignore edges going in the wrong direction */
      map[d] = e;
  }
}

bool Digger::run()
{
  Map map;
  sortEdges(mesh, snapTag, vert, edges, map);
  APF_ITERATE(Map, map, it)
    if (tryToCollapse(it->second))
      return true;
  return false;
}

}

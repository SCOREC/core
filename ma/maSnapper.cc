/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maSnapper.h"
#include "maAdapt.h"
#include "maShape.h"
#include <apfCavityOp.h>

namespace ma {

Snapper::Snapper(Adapt* a, Tag* st, bool is)
{
  adapter = a;
  snapTag = st;
  collapse.Init(a);
  isSimple = is;
  dug = false;
}

bool Snapper::setVert(Entity* v, apf::CavityOp* o)
{
  vert = v;
  if (!o->requestLocality(&vert, 1))
    return false;
  if (isSimple)
    return true;
/* in order to try an edge collapse (we don't yet know
   which edge), bring in a cavity such that all adjacent
   edges have both vertices local.
   This is basically two layers of elements around the vertex */
  apf::Up edges;
  adapter->mesh->getUp(vert,edges);
  apf::Up ovs;
  ovs.n = edges.n;
  for (int i = 0; i < edges.n; ++i)
    ovs.e[i] = apf::getEdgeVertOppositeVert(adapter->mesh, edges.e[i], vert);
  return o->requestLocality(&ovs.e[0], ovs.n);
}

static bool trySnapping(Adapt* adapter, Tag* tag, Entity* vert,
    Entity*& worstElement)
{
  Mesh* mesh = adapter->mesh;
  Vector x = getPosition(mesh, vert);
  Vector s;
  mesh->getDoubleTag(vert, tag, &s[0]);
/* move the vertex to the desired point */
  mesh->setPoint(vert, 0, s);
  Upward elements;
  mesh->getAdjacent(vert, mesh->getDimension(), elements);
  Entity* worst;
/* check resulting cavity */
  double quality = getWorstQuality(adapter, elements, &worst);
  if (quality < adapter->input->validQuality) {
    /* not ok, put the vertex back where it was */
    mesh->setPoint(vert, 0, x);
    worstElement = worst;
    return false;
  } else {
    /* ok, take off the snap tag */
    mesh->removeTag(vert, tag);
    return true;
  }
}

/* when simple snapping fails, we identify the worst element
   (most inside-out) and try to push through it by collapsing
   it. This is done by choosing one edge on that element from
   the snapping vertex and collapsing the opposite vertex
   into the snapping vertex.
   The edge is chosen to be the one most aligned with the
   desired snap vector */

static Entity* pickEdgeToDig(
    Mesh* mesh,
    Tag* snapTag,
    Entity* vert,
    Entity* worstElement)
{
  Vector x = getPosition(mesh, vert);
  Vector s;
  mesh->getDoubleTag(vert, snapTag, &s[0]);
  Vector snapVector = s - x;
  Downward es;
  int ne = mesh->getDownward(worstElement, 1, es);
  double bestDot = 0;
  Entity* edge = 0;
  for (int i = 0; i < ne; ++i) {
    Entity* vs[2];
    mesh->getDownward(es[i], 0, vs);
    if (vs[0] != vert)
      std::swap(vs[0], vs[1]);
    if (vs[0] != vert)
      continue;
    Vector ox = getPosition(mesh, vs[1]);
    Vector edgeVector = ox - x;
    double dot = edgeVector * snapVector;
    if ((!edge) || (dot > bestDot)) {
      bestDot = dot;
      edge = es[i];
    }
  }
  return edge;
}

static bool tryDigging(Adapt* adapter, Collapse& collapse,
    Entity* e, Entity* vert)
{
  Mesh* mesh = adapter->mesh;
  Entity* ov = apf::getEdgeVertOppositeVert(mesh, e, vert);
  if (!setupCollapse(collapse, e, ov))
    return false;
  double q = adapter->input->validQuality;
  if ( ! collapse.tryThisDirection(q))
    return false;
  collapse.destroyOldElements();
  return true;
}

bool Snapper::run()
{
  dug = false;
  Entity* worstElement;
  bool ok = trySnapping(adapter, snapTag, vert, worstElement);
  if (isSimple)
    return ok;
  if (!ok) {
    Entity* edge = pickEdgeToDig(adapter->mesh, snapTag, vert, worstElement);
    dug = tryDigging(adapter, collapse, edge, vert);
    if (!dug)
      return false;
    return trySnapping(adapter, snapTag, vert, worstElement);
  }
}

}

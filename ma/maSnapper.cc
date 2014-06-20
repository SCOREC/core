/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maSnapper.h"
#include "maAdapt.h"
#include "maShapeHandler.h"
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

static void collectBadElements(Adapt* a, Upward& es, apf::Up& bes)
{
  bes.n = 0;
  for (size_t i = 0; i < es.getSize(); ++i) {
    double quality = a->shape->getQuality(es[i]);
    if (quality < a->input->validQuality)
      bes.e[bes.n++] = es[i];
  }
  assert(bes.n < (int)(sizeof(bes.e) / sizeof(Entity*)));
}

static bool trySnapping(Adapt* adapter, Tag* tag, Entity* vert,
    apf::Up& badElements)
{
  Mesh* mesh = adapter->mesh;
  Vector x = getPosition(mesh, vert);
  Vector s;
  mesh->getDoubleTag(vert, tag, &s[0]);
/* move the vertex to the desired point */
  mesh->setPoint(vert, 0, s);
  Upward elements;
  mesh->getAdjacent(vert, mesh->getDimension(), elements);
/* check resulting cavity */
  collectBadElements(adapter, elements, badElements);
  if (badElements.n) {
    /* not ok, put the vertex back where it was */
    mesh->setPoint(vert, 0, x);
    return false;
  } else {
    /* ok, take off the snap tag */
    mesh->removeTag(vert, tag);
    return true;
  }
}

static bool tryDiggingEdge(Adapt* adapter, Collapse& collapse, Entity* e)
{
  Mesh* mesh = adapter->mesh;
  assert(mesh->getType(e) == EDGE);
  if ( ! collapse.setEdge(e))
    return false;
  if ( ! collapse.checkClass())
    return false;
  if ( ! collapse.checkTopo())
    return false;
  double q = adapter->input->validQuality;
  if ( ! collapse.tryBothDirections(q))
    return false;
  collapse.destroyOldElements();
  return true;
}

static bool tryDigging2(Adapt* a, Collapse& c, apf::Up& badElements)
{
  Mesh* m = a->mesh;
  for (int i = 0; i < badElements.n; ++i) {
    Entity* elem = badElements.e[i];
    Vector center = apf::getLinearCentroid(m, elem);
    Downward edges;
    int nedges = m->getDownward(elem, 1, edges);
    for (int j = 0; j < nedges; ++j)
      if (tryDiggingEdge(a, c, edges[j]))
        return true;
  }
  return false;
}

static bool tryDigging(Adapt* a, Collapse& c, Entity* v,
    apf::Up& badElements)
{
  bool hadItBefore = getFlag(a, v, DONT_COLLAPSE);
  setFlag(a, v, DONT_COLLAPSE);
  bool ok = tryDigging2(a, c, badElements);
  if (!hadItBefore)
    clearFlag(a, v, DONT_COLLAPSE);
  return ok;
}

bool Snapper::run()
{
  dug = false;
  apf::Up badElements;
  bool ok = trySnapping(adapter, snapTag, vert, badElements);
  if (isSimple)
    return ok;
  if (ok)
    return true;
  dug = tryDigging(adapter, collapse, vert, badElements);
  if (!dug)
    return false;
  return trySnapping(adapter, snapTag, vert, badElements);
}

}

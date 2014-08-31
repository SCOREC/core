#include "maLayerCollapse.h"
#include "maCrawler.h"
#include "maShape.h"

#include <cstdio>

namespace ma {

LayerCollapse::LayerCollapse(Adapt* a_)
{
  a = a_;
  m = a->mesh;
  collapse.Init(a_);
}

bool LayerCollapse::setup_(Entity* edge)
{
  edges.clear();
  curves[0].clear();
  curves[1].clear();
  elementsToCollapse.clear();
  elementsToKeep.clear();
  newSimplices.setSize(0);
  newLayer.clear();
  bool ok = collapse.setEdge(edge);
  assert(ok);
  collapse.setVerts();
  Entity* v[2] = {collapse.vertToCollapse,
                  collapse.vertToKeep};
  HasFlag p(a, COLLAPSE);
  while (true) {
    if (m->isShared(v[0]) || m->isShared(v[1]))
      return false;
    edges.push_back(edge);
    curves[0].push_back(v[0]);
    curves[1].push_back(v[1]);
    edge = getOtherEdge(m, edge, p);
    if (!edge)
      break;
    setFlag(a, edge, COLLAPSE);
    m->getDownward(edge, 0, v);
    Entity* av[2] = {v[0], curves[0].back()};
    if ( ! findUpward(m, EDGE, av)) {
      std::swap(v[0], v[1]);
      av[0] = v[0];
      assert(findUpward(m, EDGE, av));
    }
    setFlag(a, v[0], COLLAPSE);
  }
  return true;
}

bool LayerCollapse::setup(Entity* edge)
{
  if (setup_(edge))
    return true;
  unmark();
  return false;
}

void LayerCollapse::computeElementSets()
{
  int dim = m->getDimension();
  for (size_t i = 0; i < edges.size(); ++i) {
    EntityArray elements;
    m->getAdjacent(edges[i], dim, elements);
    for (size_t j = 0; j < elements.getSize(); ++j)
      elementsToCollapse.insert(elements[j]);
  }
  for (size_t i = 0; i < curves[0].size(); ++i) {
    EntityArray elements;
    m->getAdjacent(curves[0][i], dim, elements);
    for (size_t j = 0; j < elements.getSize(); ++j) {
      if ( ! elementsToCollapse.count(elements[j]))
        elementsToKeep.insert(elements[j]);
    }
  }
}

bool LayerCollapse::involvesPyramids()
{
  APF_ITERATE(EntitySet, elementsToCollapse, it)
    if (m->getType(*it) == PYRAMID)
      return true;
  APF_ITERATE(EntitySet, elementsToKeep, it)
    if (m->getType(*it) == PYRAMID)
      return true;
  return false;
}

bool LayerCollapse::checkIndividualCollapses()
{
  for (size_t i = 0; i < edges.size(); ++i) {
    assert(getFlag(a, edges[i], COLLAPSE));
    bool ok = collapse.setEdge(edges[i]);
    assert( ok );
    collapse.vertToCollapse = curves[0][i];
    collapse.vertToKeep = curves[1][i];
    if ( ! collapse.checkTopo())
      return false;
    if ( ! getFlag(a, curves[0][i], COLLAPSE))
      return false;
  }
  return true;
}

/* mostly copy-pasted from ma::rebuildElement in maMesh.cc,
   this version maps from one vector of vertices to another */
static Entity* rebuildLayerElement(
    Mesh* m,
    Entity* original,
    EntityVector& oldVerts,
    EntityVector& newVerts,
    apf::BuildCallback* cb)
{
  int type = m->getType(original);
  if (type==VERT)
  {
    for (size_t i = 0; i < oldVerts.size(); ++i)
      if (oldVerts[i] == original)
        return newVerts[i];
    return original;
  }
  int d = Mesh::typeDimension[type];
  Downward down;
  int nd = m->getDownward(original,d-1,down);
  for (int i=0; i < nd; ++i)
    down[i] = rebuildLayerElement(m, down[i], oldVerts, newVerts, cb);
  return makeOrFind(m, m->toModel(original), type, down, cb);
}

/* note - there is no solution transfer or high-order shape
   transfer in layer collapsing right now.
   that would take support in maMap.h in the form of
   not just an insideness measure but also an inverse map.
   or a user-specified transfer algorithm... */

void LayerCollapse::rebuildElements()
{
  assert(elementsToKeep.size());
  newSimplices.setSize(elementsToKeep.size());
  int nsi = 0;
  APF_ITERATE(EntitySet, elementsToKeep, it) {
    int type = m->getType(*it);
    if (apf::isSimplex(type))
      newSimplices[nsi++] = rebuildElement(
          a, *it, curves[0].back(), curves[1].back());
    else
      newLayer.push_back(rebuildLayerElement(
          m, *it, curves[0], curves[1], a->buildCallback));
  }
  newSimplices.setSize(nsi);
}

bool LayerCollapse::checkValidity(double qualityToBeat)
{
  /* it is legal in some cases for a layer collapse to
     destroy all adjacent simplices (this is the
     "no material" case when collapsing a surface edge,
     but since there are layer elements on top it is
     safe and requires no re-classification */
  if (newSimplices.getSize()) {
    double quality = getWorstQuality(a, newSimplices);
    if (quality < qualityToBeat)
      return false;
  }
  for (size_t i = 0; i < newLayer.size(); ++i)
    if ( ! isLayerElementOk(m, newLayer[i]))
      return false;
  return true;
}

void LayerCollapse::destroyOldElements()
{
  APF_ITERATE(EntitySet, elementsToCollapse, it)
    destroyElement(a, *it);
  APF_ITERATE(EntitySet, elementsToKeep, it)
    destroyElement(a, *it);
}

bool LayerCollapse::apply_(double qualityToBeat)
{
  computeElementSets();
  if (involvesPyramids()) /* avoid this complexity for now */
    return false;
  if ( ! checkIndividualCollapses())
    return false;
  rebuildElements();
  if ( ! checkValidity(qualityToBeat))
    return false;
  destroyOldElements();
  return true;
}

void LayerCollapse::unmark()
{
  for (size_t i = 1; i < edges.size(); ++i)
    clearFlag(a, edges[i], COLLAPSE);
  for (size_t i = 1; i < edges.size(); ++i)
    clearFlag(a, curves[0][i], COLLAPSE);
  if (edges.size()) {
    collapse.setEdge(edges[0]);
    collapse.unmark();
  }
}

void LayerCollapse::cancel()
{
  for (size_t i = 0; i < newSimplices.getSize(); ++i)
    destroyElement(a, newSimplices[i]);
  for (size_t i = 0; i < newLayer.size(); ++i)
    destroyElement(a, newLayer[i]);
  unmark();
}

bool LayerCollapse::apply(double qualityToBeat)
{
  if (apply_(qualityToBeat))
    return true;
  cancel();
  return false;
}

}

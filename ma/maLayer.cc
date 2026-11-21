#include "maLayer.h"
#include "maAdapt.h"
#include "maRefine.h"
#include "maShape.h"
#include <sstream>
#include <set>
#include <algorithm>
#include <pcu_util.h>
#include <lionPrint.h>

namespace ma {

static long markLayerElements(Adapt* a)
{
  Mesh* m = a->mesh;
  Entity* e;
  int meshDimension = m->getDimension();
  long n = 0;
  Iterator* it = m->begin(meshDimension);
  while ((e = m->iterate(it)))
    if ( ! apf::isSimplex(m->getType(e))) {
      setFlagOnClosure(a, e, LAYER);
      ++n;
    }
  m->end(it);
  // set LAYER flag for user defined boundary layer elements if they exist.
  Tag* layerTag = m->findTag(a->input->userDefinedLayerTagName);
  if (layerTag){
    PCU_ALWAYS_ASSERT(m->getTagType(layerTag) == apf::Mesh::INT);
    it = m->begin(meshDimension);
    while ((e = m->iterate(it))) {
      if (!m->hasTag(e, layerTag))
      	continue;
      int tagValue;
      m->getIntTag(e, layerTag, &tagValue);
      if (tagValue) {
      	setFlagOnClosure(a, e, LAYER);
      	n++;
      }
    }
  }
  n = m->getPCU()->Add<long>(n);
  a->hasLayer = (n != 0);
  if ( ! a->hasLayer)
    return 0;
  PCU_ALWAYS_ASSERT(meshDimension == 3 || meshDimension == 2);
  for (int i=0; i < 4; ++i)
    syncFlag(a,i,LAYER);
  return n;
}

void freezeLayer(Adapt* a)
{
  if ( ! a->hasLayer)
    return;
  Mesh* m = a->mesh;
  Entity* e;
  Iterator* it = m->begin(0);
  while ((e = m->iterate(it)))
    if (getFlag(a, e, LAYER))
      setFlag(a, e, DONT_COLLAPSE | DONT_SNAP);
  m->end(it);
  it = m->begin(1);
  while ((e = m->iterate(it)))
    if (getFlag(a, e, LAYER))
      setFlag(a, e, DONT_COLLAPSE | DONT_SPLIT | DONT_SWAP);
  m->end(it);
  it = m->begin(m->getDimension());
  while ((e = m->iterate(it)))
    if (getFlag(a, e, LAYER))
      setFlag(a, e, OK_QUALITY);
  m->end(it);
}

void unfreezeLayer(Adapt* a)
{
  Mesh* m = a->mesh;
  Entity* e;
  Iterator* it = m->begin(0);
  while ((e = m->iterate(it)))
    if (getFlag(a, e, LAYER))
      clearFlag(a, e, DONT_COLLAPSE | DONT_SNAP);
  m->end(it);
  it = m->begin(1);
  while ((e = m->iterate(it)))
    if (getFlag(a, e, LAYER))
      clearFlag(a, e, DONT_COLLAPSE | DONT_SPLIT | DONT_SWAP);
  m->end(it);
  it = m->begin(m->getDimension());
  while ((e = m->iterate(it)))
    if (getFlag(a, e, LAYER))
      clearFlag(a, e, OK_QUALITY);
  m->end(it);
}

void resetLayer(Adapt* a)
{
  double t0 = pcu::Time();
  long n = markLayerElements(a);
  if (!n)
    return;
  freezeLayer(a);
  double t1 = pcu::Time();
  print(a->mesh->getPCU(), "marked %ld layer elements in %f seconds", n, t1 - t0);
}

void findLayerBase(Adapt* a)
{
  Mesh* m = a->mesh;
  Iterator* it = m->begin(2);
  Entity* f;
  while ((f = m->iterate(it)))
  {
    if ((m->getType(f) == apf::Mesh::TRIANGLE) &&
        (isOnModelFace(m, f)) &&
        (m->countUpward(f) == 1) &&
        (m->getType(m->getUpward(f, 0)) == apf::Mesh::PRISM))
      setFlagOnClosure(a,f,LAYER_BASE);
  }
  m->end(it);
  for (int i=0; i < 2; ++i)
    syncFlag(a,i,LAYER_BASE);
}

/* in a prism, the growth edges connect vertex i of the
   bottom triangle to vertex i+3 of the top triangle
   (see apf::prism_edge_verts), so the next vertex in the
   stack is the opposite end of that edge */
static Entity* getNextStackVert(Mesh* m, Entity* v,
    std::set<Entity*> const& visited)
{
  apf::Adjacent prisms;
  m->getAdjacent(v, 3, prisms);
  for (size_t i = 0; i < prisms.getSize(); ++i) {
    if (m->getType(prisms[i]) != apf::Mesh::PRISM)
      continue;
    Entity* pv[6];
    m->getDownward(prisms[i], 0, pv);
    int j = apf::findIn(pv, 6, v);
    Entity* ov = pv[(j + 3) % 6];
    if (!visited.count(ov))
      return ov;
  }
  return 0;
}

static bool isAdjacentToSimplex(Mesh* m, Entity* v)
{
  apf::Adjacent elements;
  m->getAdjacent(v, m->getDimension(), elements);
  for (size_t i = 0; i < elements.getSize(); ++i)
    if (apf::isSimplex(m->getType(elements[i])))
      return true;
  return false;
}

std::vector<Entity*> getLayerVertStack(Mesh* m, Entity* vert)
{
  /* walk to one end of the stack */
  std::set<Entity*> visited;
  visited.insert(vert);
  Entity* end = vert;
  Entity* next;
  while ((next = getNextStackVert(m, end, visited))) {
    visited.insert(next);
    end = next;
  }
  /* then walk back to the other end collecting vertices in order */
  std::vector<Entity*> stack;
  visited.clear();
  stack.push_back(end);
  visited.insert(end);
  while ((next = getNextStackVert(m, stack.back(), visited))) {
    visited.insert(next);
    stack.push_back(next);
  }
  /* orient base -> top, top being next to the tets */
  if (isAdjacentToSimplex(m, stack.front()) &&
      !isAdjacentToSimplex(m, stack.back()))
    std::reverse(stack.begin(), stack.end());
  return stack;
}

void allowSplitCollapseOutsideLayer(Adapt* a)
{
  Mesh* m = a->mesh;
  Entity* e;
  Iterator* it = m->begin(1);
/* these were set by ma::refine(ma::Adapt*) and ma::coarsen(ma::Adapt*)
   for performance reasons,
   but should be disabled during shape correction so that splits and
   collapses can be used */
  while ((e = m->iterate(it)))
    if ( ! getFlag(a,e,LAYER))
      clearFlag(a,e,DONT_COLLAPSE | DONT_SPLIT);
  m->end(it);
}

void allowSplitInLayer(Adapt* a)
{
  if ( ! a->input->shouldRefineLayer)
    return;
  if ( ! a->hasLayer)
    return;
  Mesh* m = a->mesh;
  Entity* e;
  Iterator* it = m->begin(1);
  while ((e = m->iterate(it)))
    if (getFlag(a,e,LAYER))
      clearFlag(a,e,DONT_SPLIT);
  m->end(it);
  print(m->getPCU(), "allowing layer refinement");
}

void collectForLayerRefine(Refine* r)
{
  if ( ! r->adapt->input->shouldRefineLayer)
    return;
  if ( ! r->adapt->hasLayer)
    return;
/* enables collection of vertices created
   at quad centers, that mapping is needed
   to refine layer elements */
  r->shouldCollect[2] = true;
}

void checkLayerShape(Mesh* m, const char* key)
{
  double t0 = pcu::Time();
  Iterator* it = m->begin(m->getDimension());
  Entity* e;
  long n = 0;
  while ((e = m->iterate(it)))
    if ( ! apf::isSimplex(m->getType(e)))
      if (!isLayerElementOk(m, e)) {
        std::stringstream ss;
        ss.precision(15);
        ss << std::scientific;
        ss << key << ": ";
        int type = m->getType(e);
        ss << "layer " << apf::Mesh::typeName[type]
           << " at " << apf::getLinearCentroid(m, e)
           << " is unsafe to tetrahedronize\n";
        std::string s = ss.str();
        lion_oprint(1,"%s",s.c_str());
        fflush(stdout);
        ++n;
      }
  m->end(it);
  n = m->getPCU()->Add<long>(n);
  double t1 = pcu::Time();
  print(m->getPCU(), "%s: checked layer quality in %f seconds: %ld unsafe elements", key, t1 - t0, n);
}

}

#include "maLayer.h"
#include "maAdapt.h"
#include "maRefine.h"
#include "maShape.h"
#include <sstream>
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

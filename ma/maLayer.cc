#include <PCU.h>
#include "maLayer.h"
#include "maAdapt.h"
#include "maRefine.h"
#include "maShape.h"
#include <sstream>

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
  PCU_Add_Longs(&n, 1);
  a->hasLayer = (n != 0);
  if ( ! a->hasLayer)
    return 0;
  assert(meshDimension == 3);
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
  double t0 = PCU_Time();
  long n = markLayerElements(a);
  if (!n)
    return;
  freezeLayer(a);
  double t1 = PCU_Time();
  print("marked %ld layer elements in %f seconds", n, t1 - t0);
}

void findLayerBase(Adapt* a)
{
  Mesh* m = a->mesh;
  Iterator* it = m->begin(2);
  Entity* f;
  while ((f = m->iterate(it)))
  {
    if (( m->getType(f)==TRI )
      &&( isOnModelFace(m, f) )
      &&( m->countUpward(f)==1 )
      &&( m->getType(m->getUpward(f,0))==PRISM ))
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
  print("allowing layer refinement");
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

void checkLayerShape(Mesh* m)
{
  double t0 = PCU_Time();
  Iterator* it = m->begin(m->getDimension());
  Entity* e;
  while ((e = m->iterate(it)))
    if ( ! apf::isSimplex(m->getType(e)))
      if ( ! isLayerElementOk(m, e)) {
        std::stringstream ss;
        ss << "warning: layer element at "
          << apf::getLinearCentroid(m, e)
          << " is unsafe to tetrahedronize\n";
        std::string s = ss.str();
        fprintf(stderr,"%s",s.c_str());
      }
  m->end(it);
  double t1 = PCU_Time();
  print("checked layer quality in %f seconds",t1 - t0);
}

}

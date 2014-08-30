#include "maMesh.h"
#include "maAdapt.h"
#include "maLayer.h"
#include "maCoarsen.h"
#include "maCrawler.h"
#include "maLayerCollapse.h"
#include <PCU.h>

/* see maCoarsen.cc for the unstructured equivalent. */
namespace ma {

static void allowLayerToCollapse(Adapt* a)
{
  Mesh* m = a->mesh;
  for (int d = 0; d < m->getDimension(); ++d) {
    Iterator* it = m->begin(d);
    Entity* e;
    while ((e = m->iterate(it)))
      if (getFlag(a, e, LAYER))
        clearFlag(a, e, DONT_COLLAPSE);
    m->end(it);
  }
}

/* usually we would use markEntities for this,
   but this needs to be restricted to LAYER_BASE,
   and there is no DONT_COLLAPSE caching */
static long markBaseEdgesToCollapse(Adapt* a)
{
  Mesh* m = a->mesh;
  Iterator* it = m->begin(1);
  SizeField* sf = a->sizeField;
  long n = 0;
  Entity* e;
  while ((e = m->iterate(it)))
    if (getFlag(a, e, LAYER_BASE)) {
      if (sf->shouldCollapse(e)) {
        setFlag(a, e, COLLAPSE);
        ++n;
      }
    }
  m->end(it);
  PCU_Add_Longs(&n, 1);
  return n;
}

struct CurveLocalizer : public Crawler
{
  CurveLocalizer(Adapt* a_):
    Crawler(a_)
  {
    a = a_;
    m = a->mesh;
    plan = new apf::Migration(m);
    tag = m->createIntTag("ma_curve_dest", 1);
  }
  void handle(Entity* v, bool inCurve, int dest)
  {
    if (inCurve) {
      m->setIntTag(v, tag, &dest);
      EntityArray elements;
      m->getAdjacent(v, m->getDimension(), elements);
      for (size_t i = 0; i < elements.getSize(); ++i) {
        int elemDest = dest;
        if (plan->has(elements[i]))
          elemDest = std::max(elemDest, plan->sending(elements[i]));
        plan->send(elements[i], elemDest);
      }
    }
    setFlag(a, v, CHECKED);
  }
  void begin(Layer& first)
  {
    getDimensionBase(a, 0, first);
    for (size_t i = 0; i < first.size(); ++i)
      setFlag(a, first[i], CHECKED);
  }
  void end()
  {
    clearFlagFromDimension(a, CHECKED, 0);
    m->destroyTag(tag);
  }
  int getDest(Entity* v)
  {
    if (!m->hasTag(v, tag))
      return -1;
    int dest;
    m->getIntTag(v, tag, &dest);
    return dest;
  }
  Entity* crawl(Entity* v)
  {
    HasFlag p(a, CHECKED);
    Entity* ov = getOtherVert(m, v, p);
    if (!ov)
      return ov;
    int dest = getDest(v);
    handle(ov, dest != -1, dest);
    return ov;
  }
  void send(Entity* e, int to)
  {
    int dest = getDest(e);
    PCU_COMM_PACK(to, dest);
  }
  bool recv(Entity* e, int from)
  {
    int dest;
    PCU_COMM_UNPACK(dest);
    if (m->hasTag(e, tag))
      return false;
    handle(e, dest != -1, dest);
    return true;
  }
  Adapt* a;
  Mesh* m;
  int flag;
  apf::Migration* plan;
  Tag* tag;
};

static apf::Migration* planLayerCollapseMigration(Adapt* a, int d)
{
  CurveLocalizer cl(a);
  Mesh* m = a->mesh;
  Iterator* it = m->begin(1);
  Entity* e;
  while ((e = m->iterate(it)))
    if (m->getModelType(m->toModel(e)) == d) {
      Entity* v[2];
      m->getDownward(e, 0, v);
      cl.handle(v[0], true, PCU_Comm_Self());
      cl.handle(v[1], true, PCU_Comm_Self());
    }
  crawlLayers(&cl);
  return cl.plan;
}

static bool wouldEmptyParts(apf::Migration* plan)
{
  apf::Mesh* m = plan->getMesh();
  int self = PCU_Comm_Self();
  size_t sendingAway = 0;
  for (int i = 0; i < plan->count(); ++i)
    if (plan->sending(plan->get(i)) != self)
      ++sendingAway;
  bool wouldEmptyThisPart = (sendingAway == m->count(m->getDimension()));
  return PCU_Or(wouldEmptyThisPart);
}

static apf::Migration* migrateForLayerCollapse(Adapt* a, int d)
{
  apf::Migration* plan = planLayerCollapseMigration(a, d);
  /* before looking for a fix, lets just detect if this ever happens */
  assert( ! wouldEmptyParts(plan));
  return plan;
}

static void collapseLocalStacks(Adapt* a,
    int& skipCount,
    int& successCount,
    int& failureCount,
    int d)
{
  LayerCollapse c(a);
  Mesh* m = a->mesh;
  Input* in = a->input;
  Iterator* it = m->begin(1);
  Entity* e;
  while ((e = m->iterate(it)))
    if (getFlag(a, e, COLLAPSE) &&
        (m->getModelType(m->toModel(e)) == d)) {
      if (!c.setup(e)) {
        ++skipCount;
        continue;
      }
      if (c.apply(in->validQuality))
        ++successCount;
      else
        ++failureCount;
    }
  m->end(it);
}

/* usually we would use CavityOp for this,
   but migrating one element at a time to
   localize a stack takes way too long.
   Instead, we use our own localizing migration
   set up by a Crawler. */
static long collapseAllStacks(Adapt* a, int d)
{
  long allSuccesses = 0;
  int skipCount;
  do {
    migrateForLayerCollapse(a, d);
    skipCount = 0;
    int successCount = 0;
    int failureCount = 0;
    collapseLocalStacks(a, skipCount, successCount, failureCount, d);
    allSuccesses += successCount;
  } while (PCU_Or(skipCount));
  PCU_Add_Longs(&allSuccesses, 1);
  return allSuccesses;
}

bool coarsenLayer(Adapt* a)
{
  if ( ! a->hasLayer)
    return false;
  if ( ! a->input->shouldCoarsenLayer)
    return false;
  double t0 = MPI_Wtime();
  allowLayerToCollapse(a);
  findLayerBase(a);
  long count = markBaseEdgesToCollapse(a);
  if ( ! count)
    return false;
  assert(checkFlagConsistency(a,1,COLLAPSE));
  Mesh* m = a->mesh;
  long successCount = 0;
  /* LAYER_BASE is on surface only */
  for (int d = 1; d < m->getDimension(); ++d) {
    checkAllEdgeCollapses(a, d);
    findIndependentSet(a);
    successCount += collapseAllStacks(a, d);
  }
  double t1 = MPI_Wtime();
  print("coarsened %li layer edges in %f seconds",successCount,t1-t0);
  resetLayer(a);
  return true;
}

}

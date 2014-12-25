#include <PCU.h>
#include "maTetrahedronize.h"
#include "maCrawler.h"
#include "maAdapt.h"
#include "maRefine.h"
#include "maLayer.h"
#include <apfNumbering.h>
#include <apfShape.h>

namespace ma {

int getDiagonalFromFlag(Adapt* a, Entity* e)
{
  if (getFlag(a,e,DIAGONAL_1))
    return 0;
  if (getFlag(a,e,DIAGONAL_2))
    return 1;
  return -1;
}

int getFlagFromDiagonal(int diagonal)
{
  if (diagonal==0)
    return DIAGONAL_1;
  assert(diagonal==1);
  return DIAGONAL_2;
}

static bool getEdgeDirection(apf::GlobalNumbering* n, Entity* e)
{
  apf::Mesh* m = getMesh(n);
  Entity* v[2];
  m->getDownward(e,0,v);
  return apf::getNumber(n, apf::Node(v[0], 0))
         <
         apf::getNumber(n, apf::Node(v[1], 0));
}

static void chooseBaseDiagonals(Adapt* a)
{
  Mesh* m = a->mesh;
  apf::Numbering* local =
    apf::numberOwnedDimension(m, "layer_base_number", 0);
  apf::GlobalNumbering* global = apf::makeGlobal(local);
  apf::synchronize(global);
  Entity* e;
  Iterator * it = m->begin(1);
  while ((e = m->iterate(it)))
    if (getFlag(a,e,LAYER_BASE))
    {
      if (getEdgeDirection(global, e))
        setFlag(a,e,DIAGONAL_1);
      else
        setFlag(a,e,DIAGONAL_2);
    }
  m->end(it);
  apf::destroyGlobalNumbering(global);
}

static Entity* getOtherQuad(Adapt* a, Entity* e, Predicate& visited)
{
  Mesh* m = a->mesh;
  apf::Up up;
  m->getUp(e,up);
  for (int i=0; i < up.n; ++i)
  {
    Entity* of = up.e[i];
    if ((m->getType(of)==QUAD)&&
        ( ! visited(of)))
      return of;
  }
  return 0;
}

static int getQuadEdgeDiagonalBit(
    Entity* edge,
    Entity** quadEdges,
    int* directions)
{
  int i = findIn(quadEdges,4,edge);
  int i_bit = i & 1;
  int dir_bit = directions[i];
  return i_bit ^ dir_bit;
}

static Entity* flagQuad(Adapt* a, Entity* q, Entity* e)
{
  Mesh* m = a->mesh;
  int diagonal = getDiagonalFromFlag(a,e);
  assert(diagonal != -1);
  Entity* es[4];
  int ds[4];
  getFaceEdgesAndDirections(m,q,es,ds);
  diagonal ^= getQuadEdgeDiagonalBit(e,es,ds);
  setFlag(a, q, getFlagFromDiagonal(diagonal));
  e = getQuadEdgeOppositeEdge(m,q,e);
  /* bit flip going out is the opposite of bit flip
     going in   V   */
  diagonal ^= 1 ^ getQuadEdgeDiagonalBit(e,es,ds);
  setFlag(a, e, getFlagFromDiagonal(diagonal));
  return e;
}

struct QuadFlagger : public Crawler
{
  QuadFlagger(Adapt* a):Crawler(a) {}
  void begin(Layer& first)
  {
    getDimensionBase(adapter, 1, first);
  }
  Entity* crawl(Entity* e)
  {
    HasFlag p(adapter, DIAGONAL_1 | DIAGONAL_2);
    Entity* q = getOtherQuad(adapter, e, p);
    Entity* e2 = 0;
    if (q)
      e2 = flagQuad(adapter, q, e);
    clearFlag(adapter, e, DIAGONAL_1 | DIAGONAL_2);
    return e2;
  }
  void send(Entity* e, int to)
  {
    int diagonal = getDiagonalFromFlag(adapter, e);
    PCU_COMM_PACK(to, diagonal);
  }
  bool recv(Entity* e, int)
  {
    int diagonal;
    PCU_COMM_UNPACK(diagonal);
    if (getFlag(adapter, e, DIAGONAL_1 | DIAGONAL_2))
      return false;
    setFlag(adapter, e, getFlagFromDiagonal(diagonal));
    return true;
  }
};

static void flagQuadDiagonals(Adapt* a)
{
  QuadFlagger op(a);
  crawlLayers(&op);
}

static void prepareLayerToTets(Adapt* a)
{
  findLayerBase(a);
  chooseBaseDiagonals(a);
  flagQuadDiagonals(a);
}

static void addAllLayerElements(Refine* r)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  int quadCount = apf::countEntitiesOfType(m,QUAD);
  int prismCount = apf::countEntitiesOfType(m,PRISM);
  int pyramidCount = apf::countEntitiesOfType(m,PYRAMID);
  r->toSplit[2].setSize(quadCount);
  r->toSplit[3].setSize(prismCount + pyramidCount);
  Entity* e;
  Iterator* it = m->begin(2);
  int nf = 0;
  while ((e = m->iterate(it)))
    if (m->getType(e)==QUAD) {
      m->setIntTag(e, r->numberTag, &nf);
      r->toSplit[2][nf] = e;
      ++nf;
    }
  m->end(it);
  assert(static_cast<size_t>(nf) == r->toSplit[2].getSize());
  it = m->begin(3);
  int nr = 0;
  while ((e = m->iterate(it)))
    if ((m->getType(e)==PRISM)||
        (m->getType(e)==PYRAMID)) {
      m->setIntTag(e, r->numberTag, &nr);
      r->toSplit[3][nr] = e;
      ++nr;
    }
  m->end(it);
  assert(static_cast<size_t>(nr) == r->toSplit[3].getSize());
}

void tetrahedronizeCommon(Refine* r)
{
  resetCollection(r);
  collectForTransfer(r);
  collectForMatching(r);
  splitElements(r);
  processNewElements(r);
  destroySplitElements(r);
  cleanupAfter(r);
}

void tetrahedronize(Adapt* a)
{
  if ( ! a->input->shouldTurnLayerToTets)
    return;
  assert(a->hasLayer);
  double t0 = PCU_Time();
  prepareLayerToTets(a);
  Refine* r = a->refine;
  addAllLayerElements(r);
  tetrahedronizeCommon(r);
  double t1 = PCU_Time();
  print("boundary layer converted to tets in %f seconds",t1-t0);
}

/* like QuadFlagger, but just sets CHECKED to find the remaining
   delinquent quads */
struct QuadMarker : public Crawler
{
  QuadMarker(Adapt* a_):
    Crawler(a_)
  {
    a = a_;
    m = a->mesh;
  }
  void begin(Layer& first)
  {
    getDimensionBase(a, 1, first);
    for (size_t i = 0; i < first.size(); ++i)
      setFlag(a, first[i], CHECKED);
  }
  void end()
  {
    clearFlagFromDimension(a, CHECKED, 1);
  }
  Entity* crawl(Entity* e)
  {
    HasFlag p(a, CHECKED);
    Entity* q = getOtherQuad(a, e, p);
    if (!q)
      return 0;
    setFlag(a, q, CHECKED);
    Entity* oe = getQuadEdgeOppositeEdge(m, q, e);
    setFlag(a, oe, CHECKED);
    return oe;
  }
  void send(Entity*, int)
  {
  }
  bool recv(Entity* e, int)
  {
    if (getFlag(a, e, CHECKED))
      return false;
    setFlag(a, e, CHECKED);
    return true;
  }
  Adapt* a;
  Mesh* m;
};

static void markGoodQuads(Adapt* a)
{
  QuadMarker op(a);
  crawlLayers(&op);
  syncFlag(a, 2, CHECKED);
}

static void markBadQuads(Adapt* a)
{
  Mesh* m = a->mesh;
  Entity* e;
  Iterator* it = m->begin(2);
  while ((e = m->iterate(it)))
    if (m->getType(e) == QUAD) {
      if ( ! getFlag(a, e, CHECKED)) {
        setFlag(a, e, SPLIT);
        setFlag(a, e, DIAGONAL_1);
      }
    }
  m->end(it);
  clearFlagFromDimension(a, CHECKED, 2);
  assert(checkFlagConsistency(a, 2, SPLIT));
  assert(checkFlagConsistency(a, 2, DIAGONAL_1));
}

static long markBadPyramids(Adapt* a)
{
  Mesh* m = a->mesh;
  Entity* e;
  long n = 0;
  Iterator* it = m->begin(2);
  while ((e = m->iterate(it)))
    if (getFlag(a, e, SPLIT)) {
      apf::Up up;
      m->getUp(e, up);
      for (int i = 0; i < up.n; ++i) {
        Entity* elem = up.e[i];
        assert(m->getType(elem) == PYRAMID);
        setFlag(a, elem, SPLIT);
        ++n;
      }
    }
  m->end(it);
  PCU_Add_Longs(&n, 1);
  return n;
}

static int countEntitiesWithFlag(Adapt* a, int flag, int dim)
{
  Mesh* m = a->mesh;
  Iterator* it = m->begin(dim);
  Entity* e;
  int n = 0;
  while ((e = m->iterate(it)))
    if (getFlag(a, e, flag))
      ++n;
  m->end(it);
  return n;
}

static void addBadPyramids(Refine* r)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  for (int d = 2; d <= 3; ++d)
    r->toSplit[d].setSize(countEntitiesWithFlag(a, SPLIT, d));
  size_t n[4] = {};
  for (int d = 2; d <= 3; ++d) {
    Iterator* it = m->begin(d);
    Entity* e;
    while ((e = m->iterate(it)))
      if (getFlag(a, e, SPLIT))
        r->toSplit[d][n[d]++] = e;
    m->end(it);
    assert(r->toSplit[d].getSize() == n[d]);
  }
}

static long prepareLayerCleanup(Adapt* a)
{
  findLayerBase(a);
  markGoodQuads(a);
  markBadQuads(a);
  return markBadPyramids(a);
}

void cleanupLayer(Adapt* a)
{
  if (!a->hasLayer)
    return;
  if (!a->input->shouldCleanupLayer)
    return;
  double t0 = PCU_Time();
  long n = prepareLayerCleanup(a);
  if (!n) {
    print("no bad pyramids found");
    return;
  }
  Refine* r = a->refine;
  addBadPyramids(r);
  tetrahedronizeCommon(r);
  double t1 = PCU_Time();
  print("tetrahedronized %ld bad pyramids in %f seconds", n, t1-t0);
}

}

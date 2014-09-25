#include <PCU.h>
#include "maCrawler.h"
#include "maRefine.h"
#include "maLayer.h"

namespace ma {

static void preventQuadEdgeSplits(Adapt* a)
{
  Mesh* m = a->mesh;
  Iterator* it = m->begin(2);
  Entity* e;
  while ((e = m->iterate(it)))
    if (m->getType(e) == QUAD) {
      Entity* qe[4];
      m->getDownward(e, 1, qe);
      for (int i = 0; i < 4; ++i)
        setFlag(a, qe[i], DONT_SPLIT);
    }
  m->end(it);
/* other parts may have copies of these edges
   but no quads next to them */
  syncFlag(a, 1, DONT_SPLIT);
}

static void allowBaseToSplit(Adapt* a)
{
  Mesh* m = a->mesh;
  Iterator* it = m->begin(1);
  Entity* e;
  while ((e = m->iterate(it)))
    if (getFlag(a, e, LAYER_BASE))
      clearFlag(a, e, DONT_SPLIT);
  m->end(it);
}

struct SplitTagger : public Crawler
{
  SplitTagger(Adapt* a_):
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
  void handle(Entity* e, bool split)
  {
    setFlag(a, e, CHECKED);
    if (split) {
      clearFlag(a, e, DONT_SPLIT);
      setFlag(a, e, SPLIT);
    }
  }
  Entity* crawl(Entity* e)
  {
    HasFlag p(a, CHECKED);
    Entity* oe = getOtherEdge(m, e, p);
    if (!oe)
      return 0;
    handle(oe, getFlag(a, e, SPLIT));
    return oe;
  }
  void send(Entity* e, int to)
  {
    bool has = getFlag(a, e, SPLIT);
    PCU_COMM_PACK(to, has);
  }
  bool recv(Entity* e, int)
  {
    bool has;
    PCU_COMM_UNPACK(has);
    if (getFlag(a, e, CHECKED))
      return false;
    handle(e, has);
    return true;
  }
  Adapt* a;
  Mesh* m;
};

static void tagSplits(Adapt* a)
{
  SplitTagger op(a);
  crawlLayers(&op);
}

static bool getAmbiguousTri(Adapt* a, Entity* t, Entity** v)
{
  int code_index = matchEntityToTemplate(a, t, v);
  int canonical_code = tri_edge_codes[code_index];
  return canonical_code == 0x3;
}

static void disambiguateBaseTri(Adapt* a, Entity* t)
{
  Entity* v[3];
  bool isAmbiguous = getAmbiguousTri(a, t, v);
  if (!isAmbiguous)
    return;
  int flag = DIAGONAL_1;
  if (getDistance(a, v[0], v[1]) > getDistance(a, v[2], v[1]))
    flag = DIAGONAL_2;
  setFlag(a, t, flag);
}

static void disambiguateBaseTris(Adapt* a)
{
  Mesh* m = a->mesh;
  Iterator* it = m->begin(2);
  Entity* t;
  while ((t = m->iterate(it)))
    if (getFlag(a, t, LAYER_BASE))
      disambiguateBaseTri(a, t);
  m->end(it);
}

static bool isSameCurl(Entity** a, Entity** b, int n)
{
  int i = findIn(a, n, b[0]);
  int j = findIn(a, n, b[1]);
  return j == ((i + 1) % n);
}

static int getPrismTriCurl(Mesh* m, Entity* prism, Entity* tri)
{
  Entity* pf[5];
  m->getDownward(prism, 2, pf);
  Entity* pv[6];
  m->getDownward(prism, 0, pv);
  Entity* tv[3];
  m->getDownward(tri, 0, tv);
  Entity** ptv;
  if (tri == pf[0])
    ptv = pv;
  else { assert(tri == pf[4]);
    ptv = pv + 3;
  }
  return isSameCurl(ptv, tv, 3) ? 0 : 1;
}

static Entity* getOtherPrism(Mesh* m, Entity* tri, Predicate& visited)
{
  apf::Up up;
  m->getUp(tri, up);
  for (int i = 0; i < up.n; ++i)
    if ((m->getType(up.e[i]) == PRISM) &&
        (!visited(up.e[i])))
      return up.e[i];
  return 0;
}

static Entity* getOtherTri(Mesh* m, Entity* prism, Predicate& visited)
{
  Entity* f[5];
  m->getDownward(prism, 2, f);
  if (!visited(f[0]))
    return f[0];
  if (!visited(f[4]))
    return f[4];
  return 0;
}

struct Disambiguator : public Crawler
{
  Disambiguator(Adapt* a_):
    Crawler(a_)
  {
    a = a_;
    m = a->mesh;
  }
  void begin(Layer& first)
  {
    getDimensionBase(a, 2, first);
    for (size_t i = 0; i < first.size(); ++i)
      setFlag(a, first[i], CHECKED);
  }
  void end()
  {
    clearFlagFromDimension(a, CHECKED, 2);
    clearFlagFromDimension(a, CHECKED, 3);
  }
  Entity* crawl(Entity* t)
  {
    HasFlag pred(a, CHECKED);
    Entity* p = getOtherPrism(m, t, pred);
    if (!p)
      return 0;
    setFlag(a, p, CHECKED);
    Entity* ot = getOtherTri(m, p, pred);
    setFlag(a, ot, CHECKED);
    int diag = getDiagonalFromFlag(a, t);
    if (diag == -1)
      return ot;
    diag ^= getPrismTriCurl(m, p, t);
    diag ^= getPrismTriCurl(m, p, ot);
    setFlag(a, ot, getFlagFromDiagonal(diag));
    return ot;
  }
  void send(Entity* t, int to)
  {
    int diag = getDiagonalFromFlag(a, t);
    PCU_COMM_PACK(to, diag);
  }
  bool recv(Entity* t, int)
  {
    int diag;
    PCU_COMM_UNPACK(diag);
    if (getFlag(a, t, CHECKED))
      return false;
    setFlag(a, t, CHECKED);
    if (diag == -1)
      return true;
    setFlag(a, t, getFlagFromDiagonal(diag));
    return true;
  }
  Adapt* a;
  Mesh* m;
};

static void disambiguateLayerTris(Adapt* a)
{
  disambiguateBaseTris(a);
  Disambiguator op(a);
  crawlLayers(&op);
}

void setupLayerForSplit(Adapt* a)
{
  if (!a->hasLayer)
    return;
  if (!a->input->shouldRefineLayer)
    return;
  unfreezeLayer(a);
  if (!a->input->isUniform) {
    preventQuadEdgeSplits(a);
    findLayerBase(a);
    allowBaseToSplit(a);
  }
}

void setupRefineForLayer(Refine* r)
{
  Adapt* a = r->adapt;
  if (!a->hasLayer)
    return;
  if (!a->input->shouldRefineLayer)
    return;
  if (!a->input->isUniform) {
    tagSplits(a);
    disambiguateLayerTris(a);
  }
  collectForLayerRefine(r);
}

}

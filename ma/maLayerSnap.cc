#include <PCU.h>
#include "maCrawler.h"
#include "maLayer.h"
#include "maSnap.h"
#include "maShape.h"
#include <cassert>

namespace ma {

/* this class propagates desired positions
   from layer base vertices to all other
   vertices in the layer.
   all vertices belonging to the same curve
   will be moved by the same vector
   as the layer base vertex. */
struct SnapTagger : public Crawler
{
  SnapTagger(Adapt* a_, Tag* t_):
    Crawler(a_)
  {
    a = a_;
    m = a->mesh;
    snapTag = t_;
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
  }
  Entity* crawl(Entity* v)
  {
    HasFlag p(a, CHECKED);
    Entity* ov = getOtherVert(m, v, p);
    if (!ov)
      return 0;
    setFlag(a, ov, CHECKED);
    if (m->hasTag(v, snapTag)) {
      Vector x = getPosition(m, v);
      Vector ox = getPosition(m, ov);
      Vector s;
      m->getDoubleTag(v, snapTag, &s[0]);
      Vector os = ox + (s - x);
      m->setDoubleTag(ov, snapTag, &os[0]);
    } else {
      /* tagVertsToSnap can leave a tag on the
         side of a boundary layer along a model face,
         we have to clear those away if the base vertex
         is not snapping */
      m->removeTag(ov, snapTag);
    }
    return ov;
  }
  void send(Entity* v, int to)
  {
    bool has = m->hasTag(v, snapTag);
    PCU_COMM_PACK(to, has);
    if (has) {
      Vector s;
      m->getDoubleTag(v, snapTag, &s[0]);
      PCU_COMM_PACK(to, s);
    }
  }
  bool recv(Entity* v, int)
  {
    bool has;
    PCU_COMM_UNPACK(has);
    Vector s;
    if (has)
      PCU_COMM_UNPACK(s);
    if (getFlag(a, v, CHECKED))
      return false;
    setFlag(a, v, CHECKED);
    if (has)
      m->setDoubleTag(v, snapTag, &s[0]);
    return true;
  }
  Adapt* a;
  Mesh* m;
  Tag* snapTag;
};

static void tagLayerForSnap(Adapt* a, Tag* snapTag)
{
  SnapTagger op(a, snapTag);
  crawlLayers(&op);
}

/* this class tags each layer vertex
   with a remote-copy style pointer to the vertex
   at the base of its curve.
   the pointer is (peer, idx) where
   idx is a numbering of the owned layer
   base vertices only.
   the array of owned layer base vertices is stored
   in this object, so you need this object to
   trace back to the actual base vertex pointer */
struct BaseTopLinker : public Crawler
{
  BaseTopLinker(Adapt* a_):
    Crawler(a_)
  {
    a = a_;
    m = a->mesh;
    linkTag = m->createIntTag("ma_base_top", 2);
  }
  ~BaseTopLinker()
  {
    apf::removeTagFromDimension(m, linkTag, 0);
    m->destroyTag(linkTag);
  }
  void getLink(Entity* v, int& peer, int& idx)
  {
    int link[2];
    m->getIntTag(v, linkTag, link);
    peer = link[0];
    idx = link[1];
  }
  void setLink(Entity* v, int peer, int idx)
  {
    int link[2];
    link[0] = peer;
    link[1] = idx;
    m->setIntTag(v, linkTag, link);
  }
  bool hasLink(Entity* v)
  {
    return m->hasTag(v, linkTag);
  }
  void begin(Layer& first)
  {
    getDimensionBase(a, 0, first);
    int peer = PCU_Comm_Self();
    for (size_t i = 0; i < first.size(); ++i) {
      if (!m->isOwned(first[i]))
        continue;
      int idx = base.size();
      base.push_back(first[i]);
      setLink(first[i], peer, idx);
    }
    first = base;
    syncLayer(this, first);
  }
  Entity* crawl(Entity* v)
  {
    HasTag p(m, linkTag);
    Entity* ov = getOtherVert(m, v, p);
    if (!ov)
      return 0;
    int peer, idx;
    getLink(v, peer, idx);
    setLink(ov, peer, idx);
    return ov;
  }
  void send(Entity* v, int to)
  {
    int link[2];
    m->getIntTag(v, linkTag, link);
    PCU_COMM_PACK(to, link);
  }
  bool recv(Entity* v, int)
  {
    int link[2];
    PCU_COMM_UNPACK(link);
    if (hasLink(v))
      return false;
    m->setIntTag(v, linkTag, link);
    return true;
  }
  Entity* lookup(int idx)
  {
    return base[idx];
  }
  Adapt* a;
  Mesh* m;
  Tag* linkTag;
  Layer base;
};

/* for each layer curve whose base vertex
   has a snapTag, this class moves the
   curve vertices to
   their desired positions.
   it also stores their old positions in the snapTag,
   for possible subsequent unsnapping.
   at the same time, it removes the snapTag
   from curves whose base vertex does not
   have it. */
struct LayerSnapper : public Crawler
{
  LayerSnapper(Adapt* a_, Tag* t_):
    Crawler(a_)
  {
    a = a_;
    m = a->mesh;
    snapTag = t_;
  }
  void snap(Entity* v)
  {
    Vector s;
    Vector x;
    m->getDoubleTag(v, snapTag, &s[0]);
    m->getPoint(v, 0, x);
    m->setDoubleTag(v, snapTag, &x[0]); //save old spot for unsnapping
    m->setPoint(v, 0, s);
  }
  void handle(Entity* v, bool shouldSnap)
  {
    setFlag(a, v, CHECKED);
    if (shouldSnap) {
      snap(v);
    } else if (m->hasTag(v, snapTag)) {
      m->removeTag(v, snapTag);
    }
  }
  void begin(Layer& first)
  {
    ncurves = 0;
    getDimensionBase(a, 0, first);
    Layer owned;
    for (size_t i = 0; i < first.size(); ++i) {
      Entity* v = first[i];
      if (m->isOwned(v)) {
        bool isSnapping = m->hasTag(v, snapTag);
        handle(v, isSnapping);
        owned.push_back(v);
        if (isSnapping)
          ++ncurves;
      }
    }
    syncLayer(this, owned);
    PCU_Add_Longs(&ncurves, 1);
  }
  void end()
  {
    clearFlagFromDimension(a, CHECKED, 0);
  }
  Entity* crawl(Entity* v)
  {
    HasFlag p(a, CHECKED);
    Entity* ov = getOtherVert(m, v, p);
    if (!ov)
      return 0;
    handle(ov, m->hasTag(v, snapTag));
    return ov;
  }
  void send(Entity* v, int to)
  {
    bool has = m->hasTag(v, snapTag);
    PCU_COMM_PACK(to, has);
  }
  bool recv(Entity* v, int)
  {
    bool has;
    PCU_COMM_UNPACK(has);
    if (getFlag(a, v, CHECKED))
      return false;
    handle(v, has);
    return true;
  }
  Adapt* a;
  Mesh* m;
  Tag* snapTag;
  long ncurves;
};

static long snapAllCurves(Adapt* a, Tag* snapTag)
{
  double t0 = PCU_Time();
  LayerSnapper op(a, snapTag);
  crawlLayers(&op);
  double t1 = PCU_Time();
  print("snapped %ld curves in %f seconds", op.ncurves, t1 - t0);
  return op.ncurves;
}

static bool isElementOk(Adapt* a, Entity* e)
{
  if (apf::isSimplex(a->mesh->getType(e)))
    return measureTetQuality(a->mesh, a->sizeField, e) > 0;
  return isLayerElementOk(a->mesh, e);
}

/* for each layer curve whose base vertex
   has a snapTag, this class checks whether
   adjacent layer elements are OK in shape.
   if not, the LAYER_UNSNAP flag is set
   on a subset of the curve vertices including the top vertex. */
struct UnsnapChecker : public Crawler
{
  UnsnapChecker(Adapt* a_, Tag* t_):
    Crawler(a_)
  {
    a = a_;
    m = a->mesh;
    snapTag = t_;
    foundAnything = false;
  }
  void handle(Entity* v, bool alreadyUnsnapping)
  {
    setFlag(a, v, CHECKED);
    if (!m->hasTag(v, snapTag))
      return;
    if (alreadyUnsnapping) {
      setFlag(a, v, LAYER_UNSNAP);
      assert(m->hasTag(v, snapTag));
      return;
    }
    apf::Adjacent elements;
    m->getAdjacent(v, m->getDimension(), elements);
    APF_ITERATE(apf::Adjacent, elements, eit)
      if (!isElementOk(a, *eit)) {
        foundAnything = true;
        setFlag(a, v, LAYER_UNSNAP);
        assert(m->hasTag(v, snapTag));
        return;
      }
  }
  void begin(Layer& first)
  {
    getDimensionBase(a, 0, first);
    Layer owned;
    for (size_t i = 0; i < first.size(); ++i) {
      Entity* v = first[i];
      if (m->isOwned(v)) {
        handle(v, false);
        owned.push_back(v);
      }
    }
    /* see comment of crawlLayers_doubleSync */
    syncLayer(this, owned);
    syncLayer(this, owned);
  }
  void end()
  {
    clearFlagFromDimension(a, CHECKED, 0);
  }
  Entity* crawl(Entity* v)
  {
    HasFlag p(a, CHECKED);
    Entity* ov = getOtherVert(m, v, p);
    if (!ov)
      return 0;
    handle(ov, getFlag(a, v, LAYER_UNSNAP));
    return ov;
  }
  void send(Entity* v, int to)
  {
    bool has = getFlag(a, v, LAYER_UNSNAP);
    PCU_COMM_PACK(to, has);
  }
  bool recv(Entity* v, int)
  {
    bool has;
    PCU_COMM_UNPACK(has);
    bool wasChecked = getFlag(a, v, CHECKED);
    if (wasChecked) {
      if (has)
        setFlag(a, v, LAYER_UNSNAP);
      return false;
    } else {
      handle(v, has);
      return true;
    }
  }
  Adapt* a;
  Mesh* m;
  Tag* snapTag;
  bool foundAnything;
};

/* the Unsnap checker may need to sync twice.
   consider two copies of a vertex A and B.
   if the crawler initially only knows about A,
   and A is not adjacent to bad elements, it
   will transmit to B saying "add this to the layer
   but its not unsnapping".
   if B *is* adjacent to bad elements, it needs
   to return a message to A saying "actually,
   yes you are unsnapping". */
static void crawlLayers_doubleSync(Crawler* c)
{
  Crawler::Layer layer;
  c->begin(layer);
  while (PCU_Or( ! layer.empty())) {
    crawlLayer(c, layer);
    syncLayer(c, layer);
    syncLayer(c, layer);
  }
  c->end();
}

static bool checkForUnsnap(Adapt* a, Tag* snapTag)
{
  double t0 = PCU_Time();
  UnsnapChecker op(a, snapTag);
  crawlLayers_doubleSync(&op);
  bool notOk = PCU_Or(op.foundAnything);
  double t1 = PCU_Time();
  if (notOk)
    print("checked snapped curves in %f seconds, found some to unsnap", t1 - t0);
  else
    print("checked snapped curves in %f seconds, all ok", t1 - t0);
  return notOk;
}

/* for every layer top vertex that has the LAYER_UNSNAP
   flag, this flag is fed back to its base vertex. */
static void feedbackUnsnap(Adapt* a, Tag* snapTag, BaseTopLinker& l)
{
  crawlLayers(&l);
  Mesh* m = l.m;
  long n = 0;
  Entity* v;
  PCU_Comm_Begin();
  Iterator* it = m->begin(0);
  while ((v = m->iterate(it)))
    if (getFlag(a, v, LAYER_TOP) &&
        getFlag(a, v, LAYER_UNSNAP) &&
        m->isOwned(v)) {
      int peer, link;
      l.getLink(v, peer, link);
      PCU_COMM_PACK(peer, link);
      ++n;
    }
  m->end(it);
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    int link;
    PCU_COMM_UNPACK(link);
    Entity* v = l.lookup(link);
    setFlag(a, v, LAYER_UNSNAP);
    assert(m->hasTag(v, snapTag));
  }
  n = PCU_Add_Long(n);
  print("fed back unsnap flag from %ld tops", n); 
}

/* for each layer curve whose base vertex
   has the LAYER_UNSNAP flag, this class
   unsnaps that curve.
   it also removes the snapTag from unsnapped
   vertices to prevent their consideration for
   future unsnapping. */
struct Unsnapper : public Crawler
{
  Unsnapper(Adapt* a_, Tag* t_):
    Crawler(a_)
  {
    a = a_;
    m = a->mesh;
    snapTag = t_;
    ncurves = 0;
  }
  void unsnap(Entity* v)
  {
    setFlag(a, v, LAYER_UNSNAP);
    Vector s;
    m->getDoubleTag(v, snapTag, &s[0]);
    m->setPoint(v, 0, s);
    m->removeTag(v, snapTag);
  }
  void handle(Entity* v, bool shouldUnsnap)
  {
    setFlag(a, v, CHECKED);
    if (shouldUnsnap)
      unsnap(v);
  }
  void begin(Layer& first)
  {
    ncurves = 0;
    getDimensionBase(a, 0, first);
    Layer owned;
    for (size_t i = 0; i < first.size(); ++i) {
      Entity* v = first[i];
      if (m->isOwned(v)) {
        bool isUnsnapping = getFlag(a, v, LAYER_UNSNAP);
        handle(v, isUnsnapping);
        owned.push_back(v);
        if (isUnsnapping)
          ++ncurves;
      }
    }
    PCU_Add_Longs(&ncurves, 1);
    syncLayer(this, owned);
  }
  void end()
  {
    clearFlagFromDimension(a, CHECKED, 0);
    clearFlagFromDimension(a, LAYER_UNSNAP, 0);
  }
  Entity* crawl(Entity* v)
  {
    HasFlag p(a, CHECKED);
    Entity* ov = getOtherVert(m, v, p);
    if (!ov)
      return 0;
    handle(ov, getFlag(a, v, LAYER_UNSNAP));
    return ov;
  }
  void send(Entity* v, int to)
  {
    bool has = getFlag(a, v, LAYER_UNSNAP);
    PCU_COMM_PACK(to, has);
  }
  bool recv(Entity* v, int)
  {
    bool has;
    PCU_COMM_UNPACK(has);
    if (getFlag(a, v, CHECKED))
      return false;
    handle(v, has);
    return true;
  }
  Adapt* a;
  Mesh* m;
  Tag* snapTag;
  long ncurves;
};

static long unsnapMarkedCurves(Adapt* a, Tag* snapTag)
{
  double t0 = PCU_Time();
  Unsnapper op(a, snapTag);
  crawlLayers(&op);
  double t1 = PCU_Time();
  print("unsnapped %ld curves in %f seconds", op.ncurves, t1 - t0); 
  return op.ncurves;
}

void snapLayer(Adapt* a, Tag* snapTag)
{
  if ( ! a->hasLayer)
    return;
  double t0 = PCU_Time();
  findLayerBase(a);
  tagLayerForSnap(a, snapTag);
  flagLayerTop(a);
  BaseTopLinker* l = new BaseTopLinker(a);
  crawlLayers(l);
  long nsnapped = snapAllCurves(a, snapTag);
  long nunsnapped = 0;
  while (checkForUnsnap(a, snapTag)) {
    feedbackUnsnap(a, snapTag, *l);
    nunsnapped += unsnapMarkedCurves(a, snapTag);
  }
  delete l;
  double t1 = PCU_Time();
  print("finished snapping %ld of %ld layer curves in %f seconds",
      nsnapped - nunsnapped, nsnapped, t1 - t0);
}

}

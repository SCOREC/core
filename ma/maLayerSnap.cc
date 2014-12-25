#include <PCU.h>
#include "maCrawler.h"
#include "maLayer.h"
#include "maSnap.h"

namespace ma {

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

static void feedbackTopSnap(Adapt* a, Tag* snapTag)
{
  BaseTopLinker l(a);
  crawlLayers(&l);
  Mesh* m = l.m;
  Entity* v;
  PCU_Comm_Begin();
  Iterator* it = m->begin(0);
  while ((v = m->iterate(it)))
/* all layer top vertices that failed to snap */
    if (getFlag(a, v, LAYER_TOP) &&
        l.hasLink(v) &&
        m->hasTag(v, snapTag) &&
        m->isOwned(v)) {
      int peer, link;
      l.getLink(v, peer, link);
      PCU_COMM_PACK(peer, link);
    }
  m->end(it);
  PCU_Comm_Send();
  while (PCU_Comm_Receive()) {
    int link;
    PCU_COMM_UNPACK(link);
    Entity* v = l.lookup(link);
    m->removeTag(v, snapTag);
  }
}

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
    //tops should have already snapped
    if (getFlag(a, v, LAYER_TOP))
      return;
    Vector s;
    m->getDoubleTag(v, snapTag, &s[0]);
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
    getDimensionBase(a, 0, first);
    Layer owned;
    for (size_t i = 0; i < first.size(); ++i) {
      Entity* v = first[i];
      if (m->isOwned(v)) {
        handle(v, m->hasTag(v, snapTag));
        owned.push_back(v);
      }
    }
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
};

static void snapLowerLayer(Adapt* a, Tag* snapTag)
{
  LayerSnapper op(a, snapTag);
  crawlLayers(&op);
}

static long allowTopToSnap(Adapt* a, Tag* snapTag)
{
  Mesh* m = a->mesh;
  Iterator* it = m->begin(0);
  Entity* v;
  long n = 0;
  while ((v = m->iterate(it)))
    if (getFlag(a, v, LAYER_TOP) &&
        m->hasTag(v, snapTag)) {
      clearFlag(a, v, DONT_SNAP);
      if (m->isOwned(v))
        ++n;
    }
  m->end(it);
  PCU_Add_Longs(&n, 1);
  return n;
}

static void freezeTop(Adapt* a)
{
  Mesh* m = a->mesh;
  Iterator* it = m->begin(0);
  Entity* v;
  while ((v = m->iterate(it)))
    if (getFlag(a, v, LAYER_TOP))
      setFlag(a, v, DONT_SNAP);
  m->end(it);
}

void snapLayer(Adapt* a, Tag* snapTag)
{
  if ( ! a->hasLayer)
    return;
  double t0 = PCU_Time();
  findLayerBase(a);
  tagLayerForSnap(a, snapTag);
  flagLayerTop(a);
  long targets = allowTopToSnap(a, snapTag);
  long success = snapTaggedVerts(a, snapTag);
  freezeTop(a);
  feedbackTopSnap(a, snapTag);
  snapLowerLayer(a, snapTag);
  double t1 = PCU_Time();
  print("snapped %ld of %ld layer curves in %f seconds",
      success, targets, t1 - t0);
}

}

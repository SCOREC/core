#include <PCU.h>
#include "maCrawler.h"
#include "maAdapt.h"
#include "maLayer.h"
#include <apfCavityOp.h>

namespace ma {

void syncLayer(Crawler* c, Crawler::Layer& layer)
{
  Mesh* m = c->adapter->mesh;
  PCU_Comm_Begin();
  for (size_t i = 0; i < layer.size(); ++i) {
    Entity* e = layer[i];
    if (m->isShared(e)) {
      apf::Copies remotes;
      m->getRemotes(e,remotes);
      APF_ITERATE(apf::Copies,remotes,it) {
        PCU_COMM_PACK(it->first,it->second);
        c->send(e, it->first);
      }
    }
  }
  PCU_Comm_Send();
  while (PCU_Comm_Listen()) {
    int from = PCU_Comm_Sender();
    while ( ! PCU_Comm_Unpacked()) {
      Entity* e;
      PCU_COMM_UNPACK(e);
      if (c->recv(e, from))
        layer.push_back(e);
    }
  }
}

static void crawlLayer(Crawler* c, Crawler::Layer& layer)
{
  Crawler::Layer nextLayer;
  for (size_t i = 0; i < layer.size(); ++i) {
    Entity* e = layer[i];
    Entity* e2 = c->crawl(e);
    if (e2)
      nextLayer.push_back(e2);
  }
  layer.swap(nextLayer);
}

void crawlLayers(Crawler* c)
{
  Crawler::Layer layer;
  c->begin(layer);
  while (PCU_Or( ! layer.empty())) {
    crawlLayer(c, layer);
    syncLayer(c, layer);
  }
  c->end();
}

void getDimensionBase(Adapt* a, int d, Crawler::Layer& base)
{
  Mesh* m = a->mesh;
  Iterator* it = m->begin(d);
  Entity* e;
  while ((e = m->iterate(it)))
    if (getFlag(a, e, LAYER_BASE))
      base.push_back(e);
  m->end(it);
}

Entity* getOtherVert(Mesh* m, Entity* v, Predicate& visited)
{
  Upward faces;
  m->getAdjacent(v, 2, faces);
  APF_ITERATE(Upward, faces, it) {
    if (m->getType(*it) != QUAD)
      continue;
    Entity* vs[4];
    m->getDownward(*it, 0, vs);
    int i = apf::findIn(vs, 4, v);
    int j = (i + 1) % 4;
    int k = (i - 1 + 4) % 4;
    if (!visited(vs[j]))
      return vs[j];
    if (!visited(vs[k]))
      return vs[k];
  }
  return 0;
}

Entity* getOtherEdge(Mesh* m, Entity* e, Predicate& visited)
{
  Upward faces;
  m->getAdjacent(e, 2, faces);
  APF_ITERATE(Upward, faces, it) {
    if (m->getType(*it) != QUAD)
      continue;
    Entity* es[4];
    m->getDownward(*it, 1, es);
    int i = apf::findIn(es, 4, e);
    int j = (i + 2) % 4;
    if (!visited(es[j]))
      return es[j];
  }
  return 0;
}

struct Tagger
{
  void init(Mesh* m_, Tag* t_)
  {
    m = m_;
    tag = t_;
  }
  int getNumber(Entity* v)
  {
    int n;
    m->getIntTag(v, tag, &n);
    return n;
  }
  void setNumber(Entity* v, int n)
  {
    m->setIntTag(v, tag, &n);
  }
  bool hasNumber(Entity* v)
  {
    return m->hasTag(v, tag);
  }
  Mesh* m;
  Tag* tag;
};

struct LayerNumberer : public Crawler
{
  LayerNumberer(Adapt* a_):
    Crawler(a_)
  {
    a = a_;
    m = a->mesh;
    tag = m->createIntTag("ma_layer", 1);
    t.init(m, tag);
  }
  void begin(Layer& first)
  {
    getDimensionBase(a, 0, first);
    for (size_t i = 0; i < first.size(); ++i)
      t.setNumber(first[i], 0);
  }
  Entity* crawl(Entity* v)
  {
    HasTag p(m, tag);
    Entity* ov = getOtherVert(m, v, p);
    if (!ov)
      return 0;
    t.setNumber(ov, t.getNumber(v) + 1);
    return ov;
  }
  void send(Entity* v, int to)
  {
    int n = t.getNumber(v);
    PCU_COMM_PACK(to, n);
  }
  bool recv(Entity* v, int)
  {
    int n;
    PCU_COMM_UNPACK(n);
    if (t.hasNumber(v))
      return false;
    t.setNumber(v, n);
    return true;
  }
  Adapt* a;
  Mesh* m;
  Tag* tag;
  Tagger t;
};

static Tag* numberLayer(Adapt* a)
{
  LayerNumberer op(a);
  crawlLayers(&op);
  return op.tag;
}

struct TopFlagger : public apf::CavityOp
{
  TopFlagger(Adapt* a_, Tag* t_):
    apf::CavityOp(a_->mesh)
  {
    a = a_;
    m = a->mesh;
    t.init(m, t_);
  }
  Outcome setEntity(Entity* v_)
  {
    if ((!t.hasNumber(v_)) ||
        getFlag(a, v_, CHECKED))
      return SKIP;
    if (!requestLocality(&v_, 1))
      return REQUEST;
    v = v_;
    return OK;
  }
  bool isTop()
  { /* at this point all adjacent elements are local
       and quad-adjacent vertices have layer numbers.
       a top vertex is one that not adjacent to a quad
       that has vertices with higher layer numbers.

       here we choose to determine that by first looking
       across edges to see if any higher-layer vertices
       are there, and if we do find one, then make sure
       its across a quad.
       
       when stack heights differ, there can be vertices
       across an edge with higher layer number, but that
       may be a tet or pyramid edge and if thats all then
       this vertex is still the top of its curve */
    int n = t.getNumber(v);
    apf::Up es;
    m->getUp(v, es);
    for (int i = 0; i < es.n; ++i) {
      Entity* ov = apf::getEdgeVertOppositeVert(m, es.e[i], v);
      if (!t.hasNumber(ov))
        continue;
      int on = t.getNumber(ov);
      if (on > n) {
        apf::Up fs;
        m->getUp(es.e[i], fs);
        for (int j = 0; j < fs.n; ++j)
          if (m->getType(fs.e[j]) == QUAD)
            return false;
      }
    }
    return true;
  }
  void apply()
  {
    setFlag(a, v, CHECKED);
    if (isTop())
      setFlag(a, v, LAYER_TOP);
  }
  Adapt* a;
  Mesh* m;
  Entity* v;
  Tagger t;
};

void flagLayerTop(Adapt* a)
{
  Tag* layerNumbers = numberLayer(a);
  TopFlagger op(a, layerNumbers);
  op.applyToDimension(0);
  clearFlagFromDimension(a, CHECKED, 0);
  apf::removeTagFromDimension(a->mesh, layerNumbers, 0);
  a->mesh->destroyTag(layerNumbers);
}

}

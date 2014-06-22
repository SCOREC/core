#include "maCrawler.h"
#include "maAdapt.h"
#include "maLayer.h"
#include <PCU.h>

namespace ma {

static void syncLayer(Crawler* c, Crawler::Layer& layer)
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
  std::vector<Entity*> nextLayer;
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

void QuadFlagger::begin(Layer& layer)
{
  Mesh* m = adapter->mesh;
  Iterator* it = m->begin(1);
  Entity* e;
  while ((e = m->iterate(it)))
    if (getFlag(adapter, e, LAYER_BASE))
      layer.push_back(e);
  m->end(it);
}

void QuadFlagger::end()
{
}

static Entity* getOtherQuad(Adapt* a, Entity* e)
{
  Mesh* m = a->mesh;
  apf::Up up;
  m->getUp(e,up);
  for (int i=0; i < up.n; ++i)
  {
    Entity* of = up.e[i];
    if ((m->getType(of)==QUAD)&&
        ( ! getFlag(a,of,DIAGONAL_1 | DIAGONAL_2)))
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
  setFlag(a,q,diagonalToFlag(diagonal));
  e = getQuadEdgeOppositeEdge(m,q,e);
  /* bit flip going out is the opposite of bit flip
     going in   V   */
  diagonal ^= 1 ^ getQuadEdgeDiagonalBit(e,es,ds);
  setFlag(a,e,diagonalToFlag(diagonal));
  return e;
}

Entity* QuadFlagger::crawl(Entity* e)
{
  Entity* q = getOtherQuad(adapter, e);
  Entity* e2 = 0;
  if (q)
    e2 = flagQuad(adapter, q, e);
  clearFlag(adapter, e, DIAGONAL_1 | DIAGONAL_2);
  return e2;
}

void QuadFlagger::send(Entity* e, int to)
{
  int diagonal = getDiagonalFromFlag(adapter, e);
  PCU_COMM_PACK(to, diagonal);
}

bool QuadFlagger::recv(Entity* e, int from)
{
  int diagonal;
  PCU_COMM_UNPACK(diagonal);
  if (getFlag(adapter, e, DIAGONAL_1 | DIAGONAL_2))
    return false;
  setFlag(adapter, e, diagonalToFlag(diagonal));
  return true;
}

}

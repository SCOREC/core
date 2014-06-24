#include "maTetrahedronize.h"
#include "maCrawler.h"
#include "maAdapt.h"
#include "maRefine.h"
#include "maLayer.h"
#include <apfNumbering.h>
#include <apfShape.h>
#include <PCU.h>

namespace ma {

int getDiagonalFromFlag(Adapt* a, Entity* e)
{
  if (getFlag(a,e,DIAGONAL_1))
    return 0;
  if (getFlag(a,e,DIAGONAL_2))
    return 1;
  return -1;
}

static int diagonalToFlag(int diagonal)
{
  if (diagonal==0)
    return DIAGONAL_1;
  assert(diagonal==1);
  return DIAGONAL_2;
}

static bool getEdgeDirection(apf::Numbering* n, Entity* e)
{
  apf::Mesh* m = getMesh(n);
  Entity* v[2];
  m->getDownward(e,0,v);
  int o[2];
  o[0] = m->getOwner(v[0]);
  o[1] = m->getOwner(v[1]);
  if (o[0] != o[1])
    return o[0] < o[1];
  int vn[2];
  vn[0] = apf::getNumber(n,v[0],0,0);
  vn[1] = apf::getNumber(n,v[1],0,0);
  assert(vn[0] != vn[1]);
  return vn[0] < vn[1];
}

static void chooseBaseDiagonals(Adapt* a)
{
  Mesh* m = a->mesh;
  apf::Numbering* n = apf::numberOwnedNodes(
      m,"layer_base_number",apf::getLagrange(1));
  apf::synchronize(n);
  Entity* e;
  Iterator * it = m->begin(1);
  while ((e = m->iterate(it)))
    if (getFlag(a,e,LAYER_BASE))
    {
      if (getEdgeDirection(n,e))
        setFlag(a,e,DIAGONAL_1);
      else
        setFlag(a,e,DIAGONAL_2);
    }
  m->end(it);
  apf::destroyNumbering(n);
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

struct QuadFlagger : public Crawler
{
  QuadFlagger(Adapt* a):Crawler(a) {}
  void begin(Layer& first)
  {
    getDimensionBase(adapter, 1, first);
  }
  Entity* crawl(Entity* e)
  {
    Entity* q = getOtherQuad(adapter, e);
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
  bool recv(Entity* e, int from)
  {
    int diagonal;
    PCU_COMM_UNPACK(diagonal);
    if (getFlag(adapter, e, DIAGONAL_1 | DIAGONAL_2))
      return false;
    setFlag(adapter, e, diagonalToFlag(diagonal));
    return true;
  }
};

static void flagQuadDiagonals(Adapt* a)
{
  QuadFlagger op(a);
  crawlLayers(&op);
}

static void findDelinquents(Adapt* a)
{
  Mesh* m = a->mesh;
  Iterator* faces = m->begin(2);
  Entity* f;
  bool guilty = false;
  while ((f = m->iterate(faces)))
    if ((m->getType(f) == QUAD) &&
        (getDiagonalFromFlag(a, f) == -1))
    {
      guilty = true;
      apf::Up regions;
      m->getUp(f, regions);
      int pri = 0;
      int pyr = 0;
      for (int i = 0; i < regions.n; ++i)
      {
        if (m->getType(regions.e[i])==PRISM)
          ++pri;
        else if (m->getType(regions.e[i])==PYRAMID)
          ++pyr;
        else if (m->getType(regions.e[i])==TET)
          fprintf(stderr,"tet adjacent to quad !\n");
        else if (m->getType(regions.e[i])==HEX)
          fprintf(stderr,"hex adjacent to quad !\n");
        else
          fprintf(stderr,"unbelievable!\n");
      }
      int sh = m->isShared(f);
      int md = m->getModelType(m->toModel(f));
      int mt = m->getModelTag(m->toModel(f));
      Downward v;
      m->getDownward(f, 0, v);
      Vector c = (getPosition(m, v[0]) +
                  getPosition(m, v[1]) +
                  getPosition(m, v[2]) +
                  getPosition(m, v[3])) / 4;
      fprintf(stderr,"%d missed quad pri %d pyr %d nr %d sh %d md %d mt %d at %f %f %f\n",
          PCU_Comm_Self(), pri, pyr, regions.n,
          sh, md, mt, c[0], c[1], c[2]);
/* we assume they are between two pyramids, so set them to any diagonal
   and continue */
      setFlag(a, f, DIAGONAL_1);
    }
  m->end(faces);
  if (guilty)
    apf::writeOneVtkFile("quads", m);
  PCU_Barrier();
}

static void prepareLayerToTets(Adapt* a)
{
  findLayerBase(a);
  chooseBaseDiagonals(a);
  flagQuadDiagonals(a);
  findDelinquents(a);
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
  size_t nf = 0;
  while ((e = m->iterate(it)))
    if (m->getType(e)==QUAD)
      r->toSplit[2][nf++] = e;
  m->end(it);
  assert(nf == r->toSplit[2].getSize());
  it = m->begin(3);
  size_t nr = 0;
  while ((e = m->iterate(it)))
    if ((m->getType(e)==PRISM)||
        (m->getType(e)==PYRAMID))
      r->toSplit[3][nr++] = e;
  m->end(it);
  assert(nr == r->toSplit[3].getSize());
}

void tetrahedronize(Adapt* a)
{
  if ( ! a->input->shouldTurnLayerToTets)
    return;
  assert(a->hasLayer);
  double t0 = MPI_Wtime();
  prepareLayerToTets(a);
  Refine* r = a->refine;
  addAllLayerElements(r);
  resetCollection(r);
  collectForTransfer(r);
  collectForMatching(r);
  splitElements(r);
  processNewElements(r);
  destroySplitElements(r);
  cleanupAfter(r);
  double t1 = MPI_Wtime();
  print("boundary layer converted to tets in %f seconds",t1-t0);
}

}

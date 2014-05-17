#include "maLayer.h"
#include "maAdapt.h"
#include "maRefine.h"
#include "maBalance.h"
#include <apfNumbering.h>
#include <apfCavityOp.h>
#include <apfShape.h>
#include <PCU.h>

namespace ma {

bool findLayerElements(Adapt* a)
{
  Mesh* m = a->mesh;
  Entity* e;
  int meshDimension = m->getDimension();
  Iterator* it = m->begin(meshDimension);
  int isAllSimplex = 1;
  while ((e = m->iterate(it)))
    if ( ! apf::isSimplex(m->getType(e)))
    {
      isAllSimplex = 0;
      setFlagOnClosure(a,e,LAYER);
    }
  m->end(it);
  for (int i=0; i < 4; ++i)
    syncFlag(a,i,LAYER);
  PCU_Min_Ints(&isAllSimplex,1);
  assert(isAllSimplex || (meshDimension==3));
  if (isAllSimplex)
    return false;
  a->hasLayer = true;
  print("mesh has boundary layer");
  return true;
}

void findLayerBase(Adapt* a)
{
  Mesh* m = a->mesh;
  Iterator* it = m->begin(2);
  Entity* f;
  while ((f = m->iterate(it)))
  {
    if (( m->getType(f)==TRI )
      &&( isOnModelFace(m,f) )
      &&( m->countUpward(f)==1 )
      &&( m->getType(m->getUpward(f,0))==PRISM ))
      setFlagOnClosure(a,f,LAYER_BASE);
  }
  m->end(it);
  for (int i=0; i <= 2; ++i)
    syncFlag(a,i,LAYER_BASE);
}

void initLayer(Adapt* a)
{
  a->hasLayer = false;
  if (findLayerElements(a))
    findLayerBase(a);
}

void preventChangesToLayer(Adapt* a)
{
  if ( ! a->hasLayer)
    return;
  Mesh* m = a->mesh;
  Entity* e;
  Iterator* it = m->begin(0);
/* dont collapse a BL vertex into an
   unstructured vertex across an unstructured edge */
  while ((e = m->iterate(it)))
    if (getFlag(a,e,LAYER))
      setFlag(a,e,DONT_COLLAPSE);
  m->end(it);
  it = m->begin(1);
  while ((e = m->iterate(it)))
    if (getFlag(a,e,LAYER))
      setFlag(a,e,DONT_COLLAPSE | DONT_SPLIT | DONT_SWAP);
  m->end(it);
  it = m->begin(m->getDimension());
/* don't invoke quality measures on BL elements */
  while ((e = m->iterate(it)))
    if (getFlag(a,e,LAYER))
      setFlag(a,e,OK_QUALITY);
  m->end(it);
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

int diagonalToFlag(int diagonal)
{
  if (diagonal==0)
    return DIAGONAL_1;
  assert(diagonal==1);
  return DIAGONAL_2;
}

int getDiagonalFromFlag(Adapt* a, Entity* e)
{
  if (getFlag(a,e,DIAGONAL_1))
    return 0;
  if (getFlag(a,e,DIAGONAL_2))
    return 1;
  return -1;
}

/* given a two-bit code of available choices,
   with at least one of them available,
   chooses the first one available */
int choiceToFlag(int choices)
{
  return (choices & 1) ? DIAGONAL_1 : DIAGONAL_2;
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

void chooseBaseDiagonals(Adapt* a)
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

Entity* getOtherQuad(Adapt* a, Entity* e)
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

Entity* getQuadEdgeOppositeEdge(Mesh* m, Entity* q, Entity* e)
{
  Entity* edges[4];
  m->getDownward(q,1,edges);
  int i = findIn(edges,4,e);
  i = (i+2)%4;
  return edges[i];
}

int getQuadEdgeDiagonalBit(
    Entity* edge,
    Entity** quadEdges,
    int* directions)
{
  int i = findIn(quadEdges,4,edge);
  int i_bit = i & 1;
  int dir_bit = directions[i];
  return i_bit ^ dir_bit;
}

Entity* flagQuad(Adapt* a, Entity* q, Entity* e)
{
  Mesh* m = a->mesh;
  int diagonal = getDiagonalFromFlag(a,e);
  if (diagonal == -1)
  {
    fprintf(stderr,"(flagQuad) edge on model dim %d has no flag\n",
        m->getModelType(m->toModel(e)));
    abort();
  }
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

void getBaseLayer(Adapt* a, std::vector<Entity*>& base)
{
  Mesh* m = a->mesh;
  Iterator* it = m->begin(1);
  Entity* e;
  while ((e = m->iterate(it)))
    if (getFlag(a,e,LAYER_BASE))
      base.push_back(e);
  m->end(it);
}

void flagLayer(Adapt* a, std::vector<Entity*>& layer)
{
  std::vector<Entity*> nextLayer;
  for (size_t i = 0; i < layer.size(); ++i)
  {
    Entity* e = layer[i];
    Entity* q = getOtherQuad(a,e);
    if (q)
    {
      Entity* e2 = flagQuad(a,q,e);
      nextLayer.push_back(e2);
    }
    clearFlag(a,e,DIAGONAL_1 | DIAGONAL_2);
  }
  layer.swap(nextLayer);
}

void syncLayer(Adapt* a, std::vector<Entity*>& layer)
{
  Mesh* m = a->mesh;
  PCU_Comm_Begin();
  for (size_t i = 0; i < layer.size(); ++i)
  {
    Entity* e = layer[i];
    if (m->isShared(e))
    {
      int diagonal = getDiagonalFromFlag(a,e);
      if (diagonal == -1)
      {
        fprintf(stderr,"(syncLayer) edge on model dim %d has no flag\n",
            m->getModelType(m->toModel(e)));
        abort();
      }
      apf::Copies remotes;
      m->getRemotes(e,remotes);
      APF_ITERATE(apf::Copies,remotes,it)
      {
        PCU_COMM_PACK(it->first,it->second);
        PCU_COMM_PACK(it->first,diagonal);
      }
    }
  }
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
    while ( ! PCU_Comm_Unpacked())
    {
      Entity* e;
      PCU_COMM_UNPACK(e);
      int diagonal;
      PCU_COMM_UNPACK(diagonal);
      if ( ! getFlag(a,e,DIAGONAL_1 | DIAGONAL_2))
      {
        setFlag(a,e,diagonalToFlag(diagonal));
        layer.push_back(e);
      }
    }
}

void flagQuadDiagonals(Adapt* a)
{
  std::vector<Entity*> layer;
  getBaseLayer(a,layer);
  while (parallelOr( ! layer.empty()))
  {
    flagLayer(a,layer);
    syncLayer(a,layer);
  }
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

void prepareLayerToTets(Adapt* a)
{
  chooseBaseDiagonals(a);
  flagQuadDiagonals(a);
  findDelinquents(a);
}

void addAllLayerElements(Refine* r)
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

void turnLayerToTets(Adapt* a)
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
  r->shouldCollect[2] = true;
  r->shouldCollect[3] = true;
}

void flagNewLayerEntities(Refine* r)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  if ( ! a->input->shouldRefineLayer)
    return;
  if ( ! a->hasLayer)
    return;
  EntityArray& f = r->toSplit[2];
  for (size_t i=0; i < f.getSize(); ++i)
  {
    if ( ! getFlag(a,f[i],LAYER_BASE))
      continue;
    EntityArray& ne = r->newEntities[2][i];
    for (size_t j=0; j < ne.getSize(); ++j)
      setFlagOnClosure(a,ne[j],LAYER_BASE);
  }
  for (int i=0; i <= 2; ++i)
    syncFlag(a,i,LAYER_BASE);
  EntityArray& e = r->toSplit[3];
  for (size_t i=0; i < e.getSize(); ++i)
  {
    if ( ! getFlag(a,e[i],LAYER))
      continue;
    EntityArray& ne = r->newEntities[3][i];
    for (size_t j=0; j < ne.getSize(); ++j)
    {
      int t = m->getType(ne[j]);
      if ((t==PRISM)||(t==PYRAMID))
        setFlagOnClosure(a,ne[j],LAYER);
    }
  }
  for (int i=0; i < 4; ++i)
    syncFlag(a,i,LAYER);
}

}

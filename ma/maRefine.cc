/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include <PCU.h>
#include "maRefine.h"
#include "maTemplates.h"
#include "maAdapt.h"
#include "maMesh.h"
#include "maTables.h"
#include "maMatch.h"
#include "maSolutionTransfer.h"
#include "maShapeHandler.h"
#include "maSnap.h"
#include "maLayer.h"
#include <apf.h>

namespace ma {

void addEdgePreAllocation(Refine* r, Entity* e, int counts[4])
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  ++(counts[1]);
  apf::Up fs;
  m->getUp(e,fs);
  for (int fi=0; fi < fs.n; ++fi)
  {
    Entity* f = fs.e[fi];
    if ( ! getFlag(a,f,SPLIT))
    {
      setFlag(a,f,SPLIT);
      ++(counts[2]);
      apf::Up rs;
      m->getUp(f,rs);
      for (int ri=0; ri < rs.n; ++ri)
      {
        Entity* r = rs.e[ri];
        if ( ! getFlag(a,r,SPLIT))
        {
          setFlag(a,r,SPLIT);
          ++(counts[3]);
        }
      }
    }
  }
}

void allocateRefine(Refine* r, int counts[4])
{
  for (int d=1; d <= 3; ++d)
    r->toSplit[d].setSize(counts[d]);
}

void addEdgePostAllocation(Refine* refiner, Entity* e, int counts[4])
{
  Adapt* a = refiner->adapt;
  Mesh* m = a->mesh;
  refiner->toSplit[1][counts[1]]=e;
  m->setIntTag(e,refiner->numberTag,counts+1);
  ++(counts[1]);
  apf::Up fs;
  m->getUp(e,fs);
  for (int fi=0; fi < fs.n; ++fi)
  {
    Entity* f = fs.e[fi];
    if (getFlag(a,f,SPLIT))
    {
      clearFlag(a,f,SPLIT);
      refiner->toSplit[2][counts[2]]=f;
      m->setIntTag(f,refiner->numberTag,counts+2);
      ++(counts[2]);
      apf::Up rs;
      m->getUp(f,rs);
      for (int ri=0; ri < rs.n; ++ri)
      {
        Entity* r = rs.e[ri];
        if (getFlag(a,r,SPLIT))
        {
          clearFlag(a,r,SPLIT);
          refiner->toSplit[3][counts[3]]=r;
          m->setIntTag(r,refiner->numberTag,counts+3);
          ++(counts[3]);
        }
      }
    }
  }
}

static void addAllMarkedEdges(Refine* r)
{
  Adapt* a = r->adapt;
  Entity* e;
  int n[4] = {0,0,0,0};
  Mesh* m = a->mesh;
  Iterator* it = m->begin(1);
  while ((e = m->iterate(it)))
    if (getFlag(a,e,SPLIT))
      addEdgePreAllocation(r,e,n);
  m->end(it);
  allocateRefine(r,n);
  n[1]=n[2]=n[3]=0;
  it = m->begin(1);
  while ((e = m->iterate(it)))
    if (getFlag(a,e,SPLIT))
      addEdgePostAllocation(r,e,n);
  m->end(it);
}

Refine::Refine(Adapt* a)
{
  adapt = a;
  Mesh* m = a->mesh;
  numberTag = m->createIntTag("ma_refine_number",1);
  vertPlaceTag = m->createDoubleTag("ma_refine_xi",1);
}

Refine::~Refine()
{
  Mesh* m = adapt->mesh;
  m->destroyTag(numberTag);
  m->destroyTag(vertPlaceTag);
}

Entity* makeSplitVert(Refine* r, Entity* edge)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  Model* c = m->toModel(edge);
  SizeField* sf = a->sizeField;
  SolutionTransfer* st = a->solutionTransfer;
  double place = sf->placeSplit(edge);
/* placeSplit is [0,1], edge xi is [-1,1] */
  Vector xi(place*2-1,0,0);
  apf::MeshElement* me = apf::createMeshElement(m,edge);
  Vector point;
  apf::mapLocalToGlobal(me,xi,point);
  Vector param(0,0,0); //prevents uninitialized values
  if (a->input->shouldTransferParametric)
    transferParametricOnEdgeSplit(m,edge,place,param);
  Entity* vert = buildVertex(a,c,point,param);
  m->setDoubleTag(vert,r->vertPlaceTag,&(place));
  st->onVertex(me,xi,vert);
  sf->interpolate(me,xi,vert);
  apf::destroyMeshElement(me);
  return vert;
}

Entity* buildSplitElement(
    Refine* r,
    Entity* parent,
    int type,
    Entity** verts)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  return buildElement(a,m->toModel(parent),type,verts);
}

Entity* findSplitVert(Refine* r, int dimension, int id)
{
  Mesh* m = r->adapt->mesh;
  EntityArray& a = r->newEntities[dimension][id];
  for (size_t i=0; i < a.getSize(); ++i)
    if (m->getType(a[i])==VERT)
      return a[i];
  return 0;
}

Entity* findSplitVert(Refine* r, Entity* parent)
{
  Mesh* m = r->adapt->mesh;
  int n;
  m->getIntTag(parent,r->numberTag,&n);
  int d = getDimension(m,parent);
  return findSplitVert(r,d,n);
}

Entity* findSplitVert(Refine* r, Entity* v0, Entity* v1)
{
  Entity* v[2];
  v[0] = v0; v[1] = v1;
  Entity* edge = findUpward(r->adapt->mesh,EDGE,v);
  return findSplitVert(r,edge);
}

Entity* findPlacedSplitVert(Refine* r, Entity* v0, Entity* v1, double& place)
{
  Entity* v[2];
  v[0] = v0; v[1] = v1;
  Mesh* m = r->adapt->mesh;
  Entity* edge = findUpward(m,EDGE,v);
  Entity* vert = findSplitVert(r,edge);
  m->getDoubleTag(vert,r->vertPlaceTag,&place);
  int i = getDownIndex(m,edge,v0);
/* flip xi so it goes v0 -> v1 */
  if (i!=0) place = 1-place;
  return vert;
}

static int getEdgeSplitCode(Adapt* a, Entity* e)
{
  Downward edges;
  int ne = a->mesh->getDownward(e,1,edges);
  int code = 0;
  for (int i=0; i < ne; ++i)
    if (getFlag(a,edges[i],SPLIT))
      code |= (1<<i);
  return code;
}

/* based on type and edge split code, rotates
   the entity to the standard orientation for its template,
   returning the code index for the template and the
   rotated vertices. */
int matchToTemplate(int type, Entity** vi, int code, Entity** vo)
{
  CodeMatch const* table = code_match[type];
  assert(table[code].code_index != -1);
  int rotation = table[code].rotation;
  rotateEntity(type, vi, rotation, vo);
  return table[code].code_index;
}

int matchEntityToTemplate(Adapt* a, Entity* e, Entity** vo)
{
  int code = getEdgeSplitCode(a,e);
  Mesh* m = a->mesh;
  int type = m->getType(e);
  Downward vi;
  m->getDownward(e, 0, vi);
  return matchToTemplate(type, vi, code, vo);
}

static SplitFunction* all_templates[TYPES] =
{0,//vert
 edge_templates,
 tri_templates,
 quad_templates,//quad
 tet_templates,
 0,//hex
 prism_templates,//prism
 pyramid_templates //pyramid
};

void splitElement(Refine* r, Entity* e)
{
  Adapt* a = r->adapt;
  Downward v;
  int index = matchEntityToTemplate(a,e,v);
  Mesh* m = a->mesh;
  int type = m->getType(e);
  all_templates[type][index](r,e,v);
}

static void linkNewVerts(Refine* r)
{
  if (PCU_Comm_Peers()==1)
    return;
  struct { Entity* parent; Entity* vert; } message;
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  PCU_Comm_Begin();
  for (int d=1; d < m->getDimension(); ++d)
    for (size_t i=0; i < r->newEntities[d].getSize(); ++i)
    {
      Entity* e = r->toSplit[d][i];
      if ( ! m->isShared(e))
        continue;
      message.vert = findSplitVert(r,e);
      if ( ! message.vert)
        continue;
      Remotes remotes;
      m->getRemotes(e,remotes);
      APF_ITERATE(Remotes,remotes,rp)
      {
        int to = rp->first;
        message.parent = rp->second;
        PCU_COMM_PACK(to,message);
      }
    }
  PCU_Comm_Send();
  while (PCU_Comm_Listen())
  {
    int from = PCU_Comm_Sender();
    while ( ! PCU_Comm_Unpacked())
    {
      PCU_COMM_UNPACK(message);
      Entity* v = findSplitVert(r,message.parent);
      m->addRemote(v,from,message.vert);
    }
  }
}

void resetCollection(Refine* r)
{
  r->shouldCollect[0] = false;
  /* collect edge entities by default because this is the
     mechanism for finding new vertices. */
  r->shouldCollect[1] = true;
  r->shouldCollect[2] = false;
  r->shouldCollect[3] = false;
}

void collectForMatching(Refine* r)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  if (a->input->shouldHandleMatching)
    for (int d=1; d < m->getDimension(); ++d)
      r->shouldCollect[d] = true;
}

void collectForTransfer(Refine* r)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  SolutionTransfer* st = a->solutionTransfer;
  int td = st->getTransferDimension();
  td = std::min(td, a->shape->getTransferDimension());
  for (int d = td; d <= m->getDimension(); ++d)
    r->shouldCollect[d] = true;
}

void splitElements(Refine* r)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  NewEntities cb;
  for (int d=1; d <= m->getDimension(); ++d)
  {
    bool shouldCollect = r->shouldCollect[d];
    if (shouldCollect)
    {
      r->newEntities[d].setSize(r->toSplit[d].getSize());
      setBuildCallback(a,&cb);
    }
    for (size_t i=0; i < r->toSplit[d].getSize(); ++i)
    {
      Entity* e = r->toSplit[d][i];
      if (shouldCollect)
        cb.reset();
      splitElement(r,e);
      if (shouldCollect)
        cb.retrieve(r->newEntities[d][i]);
    }
    if (shouldCollect)
      clearBuildCallback(a);
  }
}

void transferElements(Refine* r)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  SolutionTransfer* st = a->solutionTransfer;
  int td = st->getTransferDimension();
  for (int d = td; d <= m->getDimension(); ++d)
    for (size_t i=0; i < r->toSplit[d].getSize(); ++i)
      st->onRefine(r->toSplit[d][i],r->newEntities[d][i]);
  td = a->shape->getTransferDimension();
  for (int d = td; d <= m->getDimension(); ++d)
    for (size_t i=0; i < r->toSplit[d].getSize(); ++i)
      a->shape->onRefine(r->toSplit[d][i],r->newEntities[d][i]);
}

void forgetNewEntities(Refine* r)
{
  for (int d=0; d <= 3; ++d)
    r->newEntities[d].setSize(0);
}

void destroySplitElements(Refine* r)
{
  Adapt* a = r->adapt;
  Mesh* m = a->mesh;
  int D = m->getDimension();
  for (size_t i=0; i < r->toSplit[D].getSize(); ++i)
    destroyElement(a,r->toSplit[D][i]);
  for (int d=1; d <= D; ++d)
    r->toSplit[d].setSize(0);
}

void cleanSplitVerts(Refine* r)
{
  Mesh* m = r->adapt->mesh;
  Tag* tag = r->vertPlaceTag;
/* only new mid-edge vertices have the placement tag */
  for (size_t i=0; i < r->newEntities[1].getSize(); ++i)
    m->removeTag(findSplitVert(r,1,i),tag);
}

struct ShouldSplit : public Predicate
{
  ShouldSplit(Adapt* a_):a(a_) {}
  bool operator()(Entity* e)
  {
    return a->sizeField->shouldSplit(e);
  }
  Adapt* a;
};

long markEdgesToSplit(Adapt* a)
{
  ShouldSplit p(a);
  return markEntities(a, 1, p, SPLIT, DONT_SPLIT);
}

void processNewElements(Refine* r)
{
  linkNewVerts(r);
  if (PCU_Comm_Peers()>1) {
    apf::stitchMesh(r->adapt->mesh);
    r->adapt->mesh->acceptChanges();
  }
  if (r->adapt->input->shouldHandleMatching)
    matchNewElements(r);
  transferElements(r);
}

void cleanupAfter(Refine* r)
{
  cleanSplitVerts(r);
  forgetNewEntities(r);
}

bool refine(Adapt* a)
{
  double t0 = PCU_Time();
  --(a->refinesLeft);
  setupLayerForSplit(a);
  long count = markEdgesToSplit(a);
  if ( ! count) {
    freezeLayer(a);
    return false;
  }
  assert(checkFlagConsistency(a,1,SPLIT));
  Refine* r = a->refine;
  resetCollection(r);
  collectForTransfer(r);
  collectForMatching(r);
  setupRefineForLayer(r);
  addAllMarkedEdges(r);
  splitElements(r);
  processNewElements(r);
  destroySplitElements(r);
  cleanSplitVerts(r);
  forgetNewEntities(r);
  double t1 = PCU_Time();
  print("refined %li edges in %f seconds",count,t1-t0);
  resetLayer(a);
  return true;
}

}

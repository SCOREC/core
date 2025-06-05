/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maCoarsen.h"
#include "maAdapt.h"
#include "maCollapse.h"
#include "maMatchedCollapse.h"
#include "maOperator.h"
#include "maDBG.h"
#include <pcu_util.h>
#include "apfShape.h"
#include <list>

namespace ma {

class CollapseChecker : public apf::CavityOp
{
  public:
    CollapseChecker(Adapt* a, int md):
      CavityOp(a->mesh,false),
      modelDimension(md)
    {
      collapse.Init(a);
    }
    virtual Outcome setEntity(Entity* e)
    {
      Adapt* a = getAdapt();
      if (( ! getFlag(a,e,COLLAPSE))||
          (getFlag(a,e,CHECKED)))
        return SKIP;
      Mesh* m = a->mesh;
      int md = m->getModelType(m->toModel(e));
      if (md!=modelDimension)
        return SKIP;
      bool ok = collapse.setEdge(e);
      PCU_ALWAYS_ASSERT(ok);
      if ( ! collapse.requestLocality(this))
        return REQUEST;
      return OK;
    }
    virtual void apply()
    {
      Adapt* a = getAdapt();
      Entity* e = collapse.edge;
      if (collapse.checkClass())
        setFlagMatched(a,e,CHECKED);
    }
    Adapt* getAdapt() {return collapse.adapt;}
  private:
    Collapse collapse;
    int modelDimension;
};

void checkAllEdgeCollapses(Adapt* a, int modelDimension)
{
  CollapseChecker checker(a,modelDimension);
  checker.applyToDimension(1);
  clearFlagFromDimension(a,CHECKED,1);
  PCU_ALWAYS_ASSERT(checkFlagConsistency(a,1,COLLAPSE));
  PCU_ALWAYS_ASSERT(checkFlagConsistency(a,0,COLLAPSE));
}

class IndependentSetFinder : public apf::CavityOp
{
  public:
    IndependentSetFinder(Adapt* a):
      CavityOp(a->mesh),
      adapt(a)
    {
      vertex = 0;
    }
    virtual Outcome setEntity(Entity* v)
    {
      if (( ! getFlag(adapt,v,COLLAPSE))||
          (getFlag(adapt,v,CHECKED)))
        return SKIP;
      if ( ! requestLocality(&v,1))
        return REQUEST;
      vertex = v;
      return OK;
    }
    virtual void apply()
    {
      if (isRequiredForMatchedEdgeCollapse(adapt,vertex))
        setFlagMatched(adapt,vertex,CHECKED);
      else
        clearFlagMatched(adapt,vertex,COLLAPSE);
    }
  protected:
    Adapt* adapt;
    Entity* vertex;
};

void findIndependentSet(Adapt* a)
{
  IndependentSetFinder finder(a);
  finder.applyToDimension(0);
  clearFlagFromDimension(a,CHECKED,0);
  PCU_ALWAYS_ASSERT(checkFlagConsistency(a, 0, COLLAPSE));
}

class AllEdgeCollapser : public Operator
{
  public:
    AllEdgeCollapser(Adapt* a, int md):
      modelDimension(md)
    {
      collapse.Init(a);
      successCount = 0;
      if(a->input->shouldForceAdaptation)
        qualityToBeat = getAdapt()->input->validQuality;
      else
        qualityToBeat = getAdapt()->input->goodQuality;
    }
    virtual int getTargetDimension() {return 1;}
    virtual bool shouldApply(Entity* e)
    {
      Adapt* a = getAdapt();
      if ( ! getFlag(a,e,COLLAPSE))
        return false;
      Mesh* m = a->mesh;
      int md = m->getModelType(m->toModel(e));
      if (md!=modelDimension)
        return false;
      bool ok = collapse.setEdge(e);
      PCU_ALWAYS_ASSERT(ok);
      return true;
    }
    virtual bool requestLocality(apf::CavityOp* o)
    {
      return collapse.requestLocality(o);
    }
    virtual void apply()
    {
      if ( ! collapse.checkTopo())
        return;
      if ( ! collapse.tryBothDirections(qualityToBeat))
        return;
      collapse.destroyOldElements();
      ++successCount;
    }
    Adapt* getAdapt() {return collapse.adapt;}
    int successCount;
  private:
    Collapse collapse;
    int modelDimension;
    double qualityToBeat;
};

int collapseAllEdges(Adapt* a, int modelDimension)
{
  AllEdgeCollapser collapser(a,modelDimension);
  applyOperator(a,&collapser);
  return collapser.successCount;
}

class MatchedEdgeCollapser : public Operator
{
  public:
    MatchedEdgeCollapser(Adapt* a, int md):
      modelDimension(md),
      collapse(a)
    {
      successCount = 0;
    }
    virtual int getTargetDimension() {return 1;}
    virtual bool shouldApply(Entity* e)
    {
      Adapt* a = getAdapt();
      if ( ! getFlag(a,e,COLLAPSE))
        return false;
      Mesh* m = a->mesh;
      int md = m->getModelType(m->toModel(e));
      if (md!=modelDimension)
        return false;
      collapse.setEdge(e);
      return true;
    }
    virtual bool requestLocality(apf::CavityOp* o)
    {
      return collapse.requestLocality(o);
    }
    virtual void apply()
    {
      double qualityToBeat = getAdapt()->input->validQuality;
      collapse.setEdges();
      if ( ! collapse.checkTopo())
        return;
      if ( ! collapse.tryBothDirections(qualityToBeat))
        return;
      collapse.destroyOldElements();
      ++successCount;
    }
    Adapt* getAdapt() {return collapse.adapt;}
    int successCount;
  private:
    int modelDimension;
    MatchedCollapse collapse;
};

static int collapseMatchedEdges(Adapt* a, int modelDimension)
{
  MatchedEdgeCollapser collapser(a, modelDimension);
  applyOperator(a, &collapser);
  return collapser.successCount;
}

struct ShouldCollapse : public Predicate
{
  ShouldCollapse(Adapt* a_):a(a_) {}
  bool operator()(Entity* e)
  {
    return a->sizeField->shouldCollapse(e);
  }
  Adapt* a;
};

long markEdgesToCollapse(Adapt* a)
{
  ShouldCollapse p(a);
  return markEntities(a, 1, p, COLLAPSE, NEED_NOT_COLLAPSE,
                      DONT_COLLAPSE | NEED_NOT_COLLAPSE);
}

bool oldcoarsen(Adapt* a)
{
  if (!a->input->shouldCoarsen)
    return false;
  double t0 = pcu::Time();
  --(a->coarsensLeft);
  long count = markEdgesToCollapse(a);
  if ( ! count)
    return false;
  Mesh* m = a->mesh;
  int maxDimension = m->getDimension();
  PCU_ALWAYS_ASSERT(checkFlagConsistency(a,1,COLLAPSE));
  long successCount = 0;
  for (int modelDimension=1; modelDimension <= maxDimension; ++modelDimension)
  {
    checkAllEdgeCollapses(a,modelDimension);
    findIndependentSet(a);
    if (m->hasMatching())
      successCount += collapseMatchedEdges(a, modelDimension);
    else
      successCount += collapseAllEdges(a, modelDimension);
  }
  successCount = m->getPCU()->Add<long>(successCount);
  double t1 = pcu::Time();
  print(m->getPCU(), "coarsened %li edges in %f seconds", successCount,t1-t0);
  return true;
}

void printIndependentSet(Adapt* a)
{
  apf::writeVtkFiles("independentMesh", a->mesh);
  ma_dbg::dumpMeshWithFlag(a, 0, 0, CHECKED, "independentVerts", "independentVerts");
  exit(0);
}

static bool tryCollapseEdge(Adapt* a, Entity* edge, Entity* keep, Collapse& collapse)
{
  PCU_ALWAYS_ASSERT(a->mesh->getType(edge) == apf::Mesh::EDGE);
  bool alreadyFlagged = true;
  if (keep) alreadyFlagged = getFlag(a, keep, DONT_COLLAPSE);
  if (!alreadyFlagged) setFlag(a, keep, DONT_COLLAPSE);

  bool result = false;
  if (collapse.setEdge(edge) && 
      collapse.checkClass() &&
      collapse.checkTopo() &&
      collapse.tryBothDirections(a->input->goodQuality)) {
    if (collapse.edgeGrewPastMaxLength()) {
      result = false;
      collapse.cancel();
    }
    else result = true;
  }  
  if (!alreadyFlagged) clearFlag(a, keep, DONT_COLLAPSE);
  return result;
}


Entity* getShortestEdge(Adapt* a, apf::Up& edges)
{
  double minLength = 999999;
  Entity* minEdge = edges.e[0];
  for (int i=0; i < edges.n; i++) {
    if (!getFlag(a, edges.e[i], COARSEN)) continue;
    double length = a->sizeField->measure(edges.e[i]);
    if (length < minLength) {
      minLength = length;
      minEdge = edges.e[i];
    }
  }
  return minEdge;
}

void flagIndependentSet(Adapt* a, apf::Up& edges, int& checked)
{
  for (int i=0; i < edges.n; i++) {
    Entity* vertices[2];
    a->mesh->getDownward(edges.e[i],0, vertices);
    for (int i = 0; i < 2; i++) {
      setFlag(a, vertices[i], NEED_NOT_COLLAPSE);
      if (getFlag(a, vertices[i], CHECKED)){
        clearFlag(a, vertices[i], CHECKED);
        checked--;
      }
    }
  }
}

void clearListFlag(Adapt* a, std::list<Entity*> list, int flag) 
{
  auto i = list.begin();
  while (i != list.end())
    clearFlag(a, *i++, flag);
}

bool isIndependent(Adapt* a, Entity* vertex)
{
  if (getFlag(a, vertex, NEED_NOT_COLLAPSE)) return false;
  apf::Up edges;
  a->mesh->getUp(vertex, edges);
  for (int i=0; i < edges.n; i++) {
    Entity* opposite = getEdgeVertOppositeVert(a->mesh, edges.e[i], vertex);
    if (getFlag(a, opposite, NEED_NOT_COLLAPSE)) return true;
  }
  return false;
}

Entity* getClosestIndependentVert(Adapt* a, std::list<Entity*>& shortEdgeVerts, std::list<Entity*>::iterator& i, bool& independentSetStarted, const int checked)
{
  while (checked < shortEdgeVerts.size())
  {
    i = shortEdgeVerts.begin();
    while (i != shortEdgeVerts.end())
    {
      Entity* vertex = *i;
      if (getFlag(a, vertex, CHECKED)) {i++; continue;}
      if (!independentSetStarted || isIndependent(a, vertex))
        return vertex;
      i++;
    }
    clearListFlag(a, shortEdgeVerts, NEED_NOT_COLLAPSE);
    independentSetStarted = false;
  }
  return 0;
}

std::list<Entity*> getShortEdgeVerts(Adapt* a)
{
  std::list<Entity*> shortEdgeVerts;
  Iterator* it = a->mesh->begin(1);
  Entity* edge;
  while ((edge = a->mesh->iterate(it))) 
  {
    if (!a->sizeField->shouldCollapse(edge)) continue; //TODO: speedup
    setFlag(a, edge, COARSEN);
    Entity* vertices[2];
    a->mesh->getDownward(edge,0,vertices);
    for (int i = 0; i < 2; i++) {
      if (getFlag(a, vertices[i], COARSEN)) continue;
      setFlag(a, vertices[i], COARSEN);
      shortEdgeVerts.push_back(vertices[i]);
    }
  }
  // ma_dbg::dumpMeshWithFlag(a, 0, 1, COARSEN, "shortEdges", "shortEdges");
  clearListFlag(a, shortEdgeVerts, COARSEN);
  return shortEdgeVerts;
}

void assertChecked(Adapt* a, std::list<Entity*>& shortEdgeVerts, const int currChecked)
{
  int realChecked = 0;
  std::list<Entity*>::iterator i = shortEdgeVerts.begin();
  while (i != shortEdgeVerts.end())
    if (getFlag(a, *i++, CHECKED)) realChecked++;
  PCU_ALWAYS_ASSERT(realChecked == currChecked);
}

bool coarsen(Adapt* a)
{
  if (!a->input->shouldCoarsen)
    return false;
  double t0 = pcu::Time();
  std::list<Entity*> shortEdgeVerts = getShortEdgeVerts(a);
  Collapse collapse;
  collapse.Init(a);
  int success = 0;
  int checked = 0;
  bool independentSetStarted = false;
  std::list<Entity*>::iterator i = shortEdgeVerts.begin();
  while (checked < shortEdgeVerts.size())
  {
    // assertChecked(a, shortEdgeVerts, checked);
    Entity* vertex = getClosestIndependentVert(a, shortEdgeVerts, i, independentSetStarted, checked);
    if (vertex == 0) continue;
    apf::Up adjacent;
    a->mesh->getUp(vertex, adjacent);
    Entity* shortEdge = getShortestEdge(a, adjacent);
    Entity* keepVertex = getEdgeVertOppositeVert(a->mesh, shortEdge, vertex);
    if (!a->sizeField->shouldCollapse(shortEdge)) {
      i = shortEdgeVerts.erase(i);
    }
    else if (tryCollapseEdge(a, shortEdge, keepVertex, collapse)) {
      flagIndependentSet(a, adjacent, checked);
      i = shortEdgeVerts.erase(i);
      independentSetStarted = true;
      success++;
      collapse.destroyOldElements();
    }
    else {
      setFlag(a, vertex, CHECKED);
      checked++;
    }
  }
  clearListFlag(a, shortEdgeVerts, CHECKED);
  ma::clearFlagFromDimension(a, NEED_NOT_COLLAPSE, 0);
  ma::clearFlagFromDimension(a, COARSEN, 1);
  double t1 = pcu::Time();
  print(a->mesh->getPCU(), "coarsened %li edges in %f seconds", success, t1-t0);
  return true;
}


}

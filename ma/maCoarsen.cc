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
#include <vector>
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

bool coarsen(Adapt* a)
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

//Measure and edge lenght and stores the result so it doesn't have to be calculated again
static double getLength(Adapt* a, Tag* lengthTag, Entity* edge)
{
  double length = 0;
  length = a->sizeField->measure(edge);
  if (a->mesh->hasTag(edge, lengthTag))
    a->mesh->getDoubleTag(edge, lengthTag, &length);
  else {
    length = a->sizeField->measure(edge);
    a->mesh->setDoubleTag(edge, lengthTag, &length);
  }
  return length;
  
}

//Make sure that a collapse will not create an edge longer than the max
bool collapseSizeCheck(Adapt* a, Entity* vertex, Entity* edge, apf::Up& adjacent)
{
  Entity* vCollapse = getEdgeVertOppositeVert(a->mesh, edge, vertex);
  for (int i=0; i<adjacent.n; i++) {
    Entity* newEdgeVerts[2]{vertex, getEdgeVertOppositeVert(a->mesh, adjacent.e[i], vCollapse)};
    Entity* newEdge = a->mesh->createEntity(apf::Mesh::EDGE, 0, newEdgeVerts);
    double length = a->sizeField->measure(newEdge);
    destroyElement(a, newEdge);
    if (length > MAXLENGTH) return false;
  }
  return true;
}

static bool tryCollapseEdge(Adapt* a, Entity* edge, Entity* keep, Collapse& collapse, apf::Up& adjacent)
{
  PCU_ALWAYS_ASSERT(a->mesh->getType(edge) == apf::Mesh::EDGE);
  bool alreadyFlagged = true;
  if (keep) alreadyFlagged = getFlag(a, keep, DONT_COLLAPSE);
  if (!alreadyFlagged) setFlag(a, keep, DONT_COLLAPSE);

  double quality = a->input->shouldForceAdaptation ? a->input->validQuality 
                                                  : a->input->goodQuality;

  bool result = false;
  if (collapse.setEdge(edge) && 
      collapse.checkClass() &&
      collapse.checkTopo() &&
      collapseSizeCheck(a, keep, edge, adjacent) &&
      collapse.tryBothDirections(quality)) {
    result = true;
  }  
  if (!alreadyFlagged) clearFlag(a, keep, DONT_COLLAPSE);
  return result;
}

/*
  Used after collpasing a vertex to flag adjacent vertices NEED_NOT_COLLAPSE, will create a
  set of collapses were no two adjacent vertices collapsed called an independent set. Will
  also clear adjacent vertices for collapse in next independent set since they might succeed now.
*/
void flagIndependentSet(Adapt* a, apf::Up& adjacent, size_t& checked)
{
  for (int adj=0; adj < adjacent.n; adj++) {
    Entity* vertices[2];
    a->mesh->getDownward(adjacent.e[adj],0, vertices);
    for (int v = 0; v < 2; v++) {
      setFlag(a, vertices[v], NEED_NOT_COLLAPSE);
      if (getFlag(a, vertices[v], CHECKED)){
        clearFlag(a, vertices[v], CHECKED); //needs to be checked again in next independent set
        checked--;
      }
    }
  }
}

struct EdgeLength
{
  Entity* edge;
  double length;
  bool operator<(const EdgeLength& other) const {
    return length < other.length;
  }
};

/*
  Given an iterator pointing to a vertex we will collapse the shortest adjacent edge and try the next
  shorted until one succeeds and then it will expand independent set. In Li's thesis it only attempts
  to collapse the shortest edge, but this gave us better results.
*/
bool collapseShortest(Adapt* a, Collapse& collapse, std::list<Entity*>& shortEdgeVerts, std::list<Entity*>::iterator& itr, size_t& checked, apf::Up& adjacent, Tag* lengthTag)
{
  Entity* vertex = *itr;
  std::vector<EdgeLength> sorted;
  for (int i=0; i < adjacent.n; i++) {
    double length = getLength(a, lengthTag, adjacent.e[i]);
    EdgeLength measured{adjacent.e[i], length};
    if (measured.length > MINLENGTH) continue;
    auto pos = std::lower_bound(sorted.begin(), sorted.end(), measured);
    sorted.insert(pos, measured);
  }
  if (sorted.size() == 0) { //performance optimization, will rarely result in a missed edge
    itr = shortEdgeVerts.erase(itr);
    return false;
  }
  for (size_t i=0; i < sorted.size(); i++) {
    Entity* keepVertex = getEdgeVertOppositeVert(a->mesh, sorted[i].edge, vertex);
    if (!tryCollapseEdge(a, sorted[i].edge, keepVertex, collapse, adjacent)) continue;
    flagIndependentSet(a, adjacent, checked);
    itr = shortEdgeVerts.erase(itr);
    collapse.destroyOldElements();
    return true;
  }
  setFlag(a, vertex, CHECKED);
  checked++;
  return false;
}

void clearListFlag(Adapt* a, std::list<Entity*> list, int flag) 
{
  auto i = list.begin();
  while (i != list.end())
    clearFlag(a, *i++, flag);
}

/*
  Iterates through shortEdgeVerts until it finds a vertex that is adjacent to an 
  independent set. We want our collapses to touch the independent set in order to
  reduce adjacent collapses, since collapsing adjacent vertices will result in a
  lower quality mesh.
*/
bool getAdjIndependentSet(Adapt* a, std::list<Entity*>& shortEdgeVerts, std::list<Entity*>::iterator& itr, bool& independentSetStarted, apf::Up& adjacent)
{
  size_t numItr=0;
  do {
    numItr++;
    if (itr == shortEdgeVerts.end()) itr = shortEdgeVerts.begin();
    Entity* vertex = *itr;
    if (getFlag(a, vertex, CHECKED)) {itr++; continue;} //Already tried to collapse
    a->mesh->getUp(vertex, adjacent);
    if (!independentSetStarted) return true;
    if (getFlag(a, vertex, NEED_NOT_COLLAPSE)) {itr++; continue;} //Too close to last collapse
    for (int i=0; i < adjacent.n; i++)
    {
      Entity* opposite = getEdgeVertOppositeVert(a->mesh, adjacent.e[i], vertex);
      if (getFlag(a, opposite, NEED_NOT_COLLAPSE)) return true; //Touching independent set
    }
    itr++;
  } while (numItr < shortEdgeVerts.size());
  clearListFlag(a, shortEdgeVerts, NEED_NOT_COLLAPSE);
  independentSetStarted = false;
  return false;
}

//returns a list of vertices that bound a short edge and flags them
std::list<Entity*> getShortEdgeVerts(Adapt* a, Tag* lengthTag)
{
  std::list<Entity*> shortEdgeVerts;
  Iterator* it = a->mesh->begin(1);
  Entity* edge;
  while ((edge = a->mesh->iterate(it))) 
  {
    double length = getLength(a, lengthTag, edge);
    if (length > MINLENGTH) continue;
    Entity* vertices[2];
    a->mesh->getDownward(edge,0,vertices);
    for (int i = 0; i < 2; i++) {
      if (getFlag(a, vertices[i], CHECKED)) continue;
      setFlag(a, vertices[i], CHECKED);
      shortEdgeVerts.push_back(vertices[i]);
    }
  }
  clearListFlag(a, shortEdgeVerts, CHECKED);
  a->mesh->end(it);
  return shortEdgeVerts;
}

/*
  Follows the alogritm in Li's thesis in order to coarsen all short edges 
  in a mesh while maintaining a decent quality mesh.
*/
bool coarsenMultiple(Adapt* a)
{
  if (!a->input->shouldCoarsen)
    return false;
  double t0 = pcu::Time();
  Tag* lengthTag = a->mesh->createDoubleTag("edge_length", 1);
  std::list<Entity*> shortEdgeVerts = getShortEdgeVerts(a, lengthTag);

  Collapse collapse;
  collapse.Init(a);
  int success = 0;
  size_t checked = 0;
  bool independentSetStarted = false;
  std::list<Entity*>::iterator itr = shortEdgeVerts.begin();
  while (checked < shortEdgeVerts.size())
  {
    apf::Up adjacent;
    if (!getAdjIndependentSet(a, shortEdgeVerts, itr, independentSetStarted, adjacent)) continue;
    if (collapseShortest(a, collapse, shortEdgeVerts, itr, checked, adjacent, lengthTag)) {
      independentSetStarted=true;
      success++;
    }
  }
  ma::clearFlagFromDimension(a, NEED_NOT_COLLAPSE | CHECKED, 0);
  a->mesh->destroyTag(lengthTag);
  double t1 = pcu::Time();
  print(a->mesh->getPCU(), "coarsened %d edges in %f seconds", success, t1-t0);
  return true;
}


}

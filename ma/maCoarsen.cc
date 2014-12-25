/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include <PCU.h>
#include "maCoarsen.h"
#include "maAdapt.h"
#include "maCollapse.h"
#include "maOperator.h"

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
      assert(ok);
      if ( ! collapse.requestLocality(this))
        return REQUEST;
      return OK;
    }
    virtual void apply()
    {
      Adapt* a = getAdapt();
      Entity* e = collapse.edge;
      if (collapse.checkClass())
        setFlag(a,e,CHECKED);
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
}

class IndependentSetFinder : public apf::CavityOp
{
  public:
    IndependentSetFinder(Adapt* a):
      CavityOp(a->mesh),
      adapt(a)
    {}
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
      if (isRequiredForAnEdgeCollapse(adapt,vertex))
        setFlag(adapt,vertex,CHECKED);
      else
        clearFlag(adapt,vertex,COLLAPSE);
    }
  private:
    Adapt* adapt;
    Entity* vertex;
};

void findIndependentSet(Adapt* a)
{
  IndependentSetFinder finder(a);
  finder.applyToDimension(0);
  clearFlagFromDimension(a,CHECKED,0);
}

class AllEdgeCollapser : public Operator
{
  public:
    AllEdgeCollapser(Adapt* a, int md):
      modelDimension(md)
    {
      collapse.Init(a);
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
      bool ok = collapse.setEdge(e);
      assert(ok);
      return true;
    }
    virtual bool requestLocality(apf::CavityOp* o)
    {
      return collapse.requestLocality(o);
    }
    virtual void apply()
    {
      double qualityToBeat = getAdapt()->input->validQuality;
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
};

static int collapseAllEdges(Adapt* a, int modelDimension)
{
  AllEdgeCollapser collapser(a,modelDimension);
  applyOperator(a,&collapser);
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
  return markEntities(a, 1, p, COLLAPSE, DONT_COLLAPSE);
}

bool coarsen(Adapt* a)
{
  double t0 = PCU_Time();
  --(a->coarsensLeft);
  long count = markEdgesToCollapse(a);
  if ( ! count)
    return false;
  Mesh* m = a->mesh;
  int maxDimension = m->getDimension();
  assert(checkFlagConsistency(a,1,COLLAPSE));
  long successCount = 0;
  for (int modelDimension=1; modelDimension <= maxDimension; ++modelDimension)
  {
    checkAllEdgeCollapses(a,modelDimension);
    findIndependentSet(a);
    successCount += collapseAllEdges(a,modelDimension);
  }
  PCU_Add_Longs(&successCount,1);
  double t1 = PCU_Time();
  print("coarsened %li edges in %f seconds",successCount,t1-t0);
  return true;
}

}

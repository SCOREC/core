
/******************************************************************************

  Copyright 2013 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/
#include "maShape.h"
#include "maSize.h"
#include "maAdapt.h"
#include "maSnap.h"
#include "maSnapper.h"
#include "maOperator.h"
#include "maEdgeSwap.h"
#include "maDoubleSplitCollapse.h"
#include "maSingleSplitCollapse.h"
#include "maFaceSplitCollapse.h"
#include "maShortEdgeRemover.h"
#include "maShapeHandler.h"
#include "maBalance.h"
#include "maDBG.h"
#include <pcu_util.h>

namespace ma {

struct IsBadQuality : public Predicate
{
  IsBadQuality(Adapt* a_):a(a_) {}
  bool operator()(Entity* e)
  {
    return a->shape->getQuality(e) < a->input->goodQuality;
  }
  Adapt* a;
};

int markBadQuality(Adapt* a)
{
  IsBadQuality p(a);
  return markEntities(a, a->mesh->getDimension(), p, BAD_QUALITY, OK_QUALITY);
}

void unMarkBadQuality(Adapt* a)
{
  Mesh* m = a->mesh;
  Iterator* it;
  Entity* e;
  it = m->begin(m->getDimension());
  while ((e = m->iterate(it))) {
    if (getFlag(a, e, ma::BAD_QUALITY))
      clearFlag(a, e, ma::BAD_QUALITY);
  }
  m->end(it);
}

double getMinQuality(Adapt* a)
{
  PCU_ALWAYS_ASSERT(a);
  Mesh* m;
  m = a->mesh;
  PCU_ALWAYS_ASSERT(m);
  Iterator* it = m->begin(m->getDimension());
  Entity* e;
  double minqual = 1;
  while ((e = m->iterate(it))) {
    if (!apf::isSimplex(m->getType(e)))
      continue;
    double qual = a->shape->getQuality(e);
    if (qual < minqual)
      minqual = qual;
  }
  m->end(it);
  return m->getPCU()->Min<double>(minqual);
}

class LargeAngleTriFixer : public Operator
{
  public:
    LargeAngleTriFixer(Adapt* a)
    {
      adapter = a;
      mesh = a->mesh;
      edgeSwap = makeEdgeSwap(a);
      ns = nf = 0;
      tri = 0;
      edge = 0;
    }
    virtual ~LargeAngleTriFixer()
    {
      delete edgeSwap;
    }
    virtual int getTargetDimension() {return 2;}
    virtual bool shouldApply(Entity* e)
    {
      if ( ! getFlag(adapter,e,BAD_QUALITY))
        return false;
      tri = e;
      // get the metric Q for angle computations
      SizeField* sf = adapter->sizeField;
      Matrix Q;
      apf::MeshElement* me = apf::createMeshElement(mesh, tri);
      Vector center(1./3.,1./3.,1./3.);
      sf->getTransform(me,center,Q);
      apf::destroyMeshElement(me);

      // pick the edge opposite to the largest angle (in metric) for swap
      Entity* edges[3];
      mesh->getDownward(e,1,edges);
      double minCos = 1.0;
      for (int i = 0; i < 3; i++) {
        Entity* current = edges[i%3];
        Entity*    next = edges[(i+1)%3];
        double cosAngle = apf::computeCosAngle(mesh, tri, current, next, Q);
        if (cosAngle < minCos) {
          minCos = cosAngle;
          edge = edges[(i+2)%3];
	}
      }
      return true;
    }
    virtual bool requestLocality(apf::CavityOp* o)
    {
      return o->requestLocality(&edge,1);
    }
    virtual void apply()
    {
        if (edgeSwap->run(edge))
        {
          ++ns;
          return;
        }
      ++nf;
      clearFlag(adapter,tri,BAD_QUALITY);
    }
  private:
    Adapt* adapter;
    Mesh* mesh;
    Entity* tri;
    Entity* edge;
    EdgeSwap* edgeSwap;
    int ns;
    int nf;
};

class QualityImprover2D : public Operator
{
  public:
    QualityImprover2D(Adapt* a)
    {
      adapter = a;
      mesh = a->mesh;
      edgeSwap = makeEdgeSwap(a);
      ns = nf = 0;
      edge = 0;
    }
    virtual ~QualityImprover2D()
    {
      delete edgeSwap;
    }
    virtual int getTargetDimension() {return 1;}
    virtual bool shouldApply(Entity* e)
    {
      if ( getFlag(adapter,e,DONT_SWAP))
        return false;
      if ( mesh->isShared(e) )
      	return false;
      edge = e;
      return true;
    }
    virtual bool requestLocality(apf::CavityOp* o)
    {
      return o->requestLocality(&edge,1);
    }
    virtual void apply()
    {
      if (edgeSwap->run(edge))
      {
	++ns;
	return;
      }
      ++nf;
    }
  private:
    Adapt* adapter;
    Mesh* mesh;
    Entity* edge;
    EdgeSwap* edgeSwap;
    int ns;
    int nf;
};

static void fixLargeAngleTets(Adapt* a)
{
  FaceSplitCollapse faceSplitCollapse(a);
  SingleSplitCollapse splitCollapse(a);
  DoubleSplitCollapse doubleSplitCollapse(a);
  EdgeSwap* edgeSwap = makeEdgeSwap(a);

  Entity* tet;
  Iterator* it = a->mesh->begin(3);
  while ((tet = a->mesh->iterate(it))) {
    if (!getFlag(a, tet, BAD_QUALITY)) continue;
    clearFlag(a, tet, BAD_QUALITY);
    Entity* verts[4];
    a->mesh->getDownward(tet, 0, verts);
    Entity* vert = verts[0];
    Entity* face = getTetFaceOppositeVert(a->mesh, tet, vert);
    Entity* ents[4];
    double area[4];
    int bit = getTetStats(a, vert, face, tet, ents, area);

    // two large dihedral angles -> key problem: two mesh edges
    if (bit==3 || bit==5 || bit==6) {
      if (edgeSwap->run(ents[0])) continue;
      if (edgeSwap->run(ents[1])) continue;
      if (splitCollapse.run(ents[0], vert)) continue;
      if (splitCollapse.run(ents[1], vert)) continue;
      if (doubleSplitCollapse.run(ents)) continue;
    }
    // three large dihedral angles -> key entity: a mesh face
    else {
      Entity* edges[3];
      a->mesh->getDownward(ents[0], 1, edges);
      if (edgeSwap->run(edges[0])) continue;
      if (edgeSwap->run(edges[1])) continue;
      if (edgeSwap->run(edges[2])) continue;
      //TODO: RUN FACE SWAP HERE
      if (faceSplitCollapse.run(ents[0], tet)) continue;
    }
  }
}

static void fixLargeAngleTris(Adapt* a)
{
  LargeAngleTriFixer fixer(a);
  applyOperator(a,&fixer);
}

static void improveQualities2D(Adapt* a)
{
  QualityImprover2D improver(a);
  applyOperator(a, &improver);
}

static double fixLargeAngles(Adapt* a)
{
  double t0 = pcu::Time();
  if (a->mesh->getDimension()==3)
    fixLargeAngleTets(a);
  else
    fixLargeAngleTris(a);
  double t1 = pcu::Time();
  return t1 - t0;
}

double improveQualities(Adapt* a)
{
  double t0 = pcu::Time();
  if (a->mesh->getDimension() == 3)
    return 0; // TODO: implement this for 3D
  else
    improveQualities2D(a);
  double t1 = pcu::Time();
  return t1 - t0;
}

void fixElementShapes(Adapt* a)
{
  if ( ! a->input->shouldFixShape)
    return;
  double t0 = pcu::Time();
  int count = markBadQuality(a);
  print(a->mesh->getPCU(), "--iter %d of shape correction loop: #bad elements %d", 0, count);
  int originalCount = count;
  int prev_count;
  double time;
  int iter = 0;
  do {
    if ( ! count)
      break;
    prev_count = count;
    time = fixLargeAngles(a);
    /* We need to snap the new verts as soon as they are
     * created (to avoid future problems). At the moment
     * new verts are created only during 3D mesh adapt, so
     * we only run a bulk snap operation if the mesh is 3D.
     */
    if (a->mesh->getDimension() == 3)
      snap(a);
    count = markBadQuality(a);
    if (count >= prev_count)
      unMarkBadQuality(a); // to make sure markEntities does not complain!
    // balance the mesh to avoid empty parts
    midBalance(a);
    iter++;
  } while(count < prev_count);
  double t1 = pcu::Time();
  print(a->mesh->getPCU(), "bad shapes down from %d to %d in %f seconds", 
        originalCount,count,t1-t0);
}

void printQuality(Adapt* a)
{
  if ( ! a->input->shouldPrintQuality)
    return;
  double minqual = getMinQuality(a);
  print(a->mesh->getPCU(), "worst element quality is %e",  minqual);
}

}

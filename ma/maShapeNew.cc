
/******************************************************************************

  Copyright 2013 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/
#include "maShapeNew.h"
#include "maSize.h"
#include "maAdapt.h"
#include "maSnap.h"
#include "maSnapper.h"
#include "maOperator.h"
#include "maEdgeSwap.h"
#include "maDoubleSplitCollapse.h"
#include "maSingleSplitCollapse.h"
#include "maFaceSplitCollapse.h"
#include "maFaceSwap.h"
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

int markBadQualityNew(Adapt* a)
{
  IsBadQuality p(a);
  return markEntities(a, a->mesh->getDimension(), p, BAD_QUALITY, OK_QUALITY);
}

void unMarkBadQualityNew(Adapt* a)
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

class FixShape
{
  public:
  Adapt* a;
  Collapse collapse;
  SingleSplitCollapse splitCollapse;
  DoubleSplitCollapse doubleSplitCollapse;
  FaceSplitCollapse faceSplitCollapse;
  EdgeSwap* edgeSwap;

  FixShape(Adapt* adapt) : splitCollapse(adapt), doubleSplitCollapse(adapt), faceSplitCollapse(adapt)
  {
    a = adapt;
    collapse.Init(a);
    edgeSwap = makeEdgeSwap(a);
  }

  bool collapseEdge(Entity* edge)
  {
    if (collapse.setEdge(edge) && 
          collapse.checkClass() &&
          collapse.checkTopo() &&
          collapse.tryBothDirections(a->input->validQuality)) {
      collapse.destroyOldElements();
      return true;
    }
    return false;
  }

  bool shortEdgeCase(Entity* tet)
  {
    Entity* edges[6];
    a->mesh->getDownward(tet, 1, edges);
    for (int i=0; i<6; i++)
      if (a->sizeField->measure(edges[i]) < MINLENGTH)
        if (collapseEdge(edges[i]))
          return true;
    return false;
  }

  double getWorstTriangle(Entity* tet, Entity*& worstTriangle)
  {
    Entity* triangles[4];
    a->mesh->getDownward(tet, 2, triangles);
    double worstQuality = a->sizeField->measure(triangles[0]);
    worstTriangle = triangles[0];
    for (int i=1; i<4; i++) {
      double quality = a->sizeField->measure(triangles[i]);
      if (quality < worstQuality) {
        worstQuality = quality;
        worstTriangle = triangles[i];
      }
    }
    return worstQuality;
  }

  bool oneLargeAngle(Entity* tet)
  {
    Entity* worstTriangle;
    if (getWorstTriangle(tet, worstTriangle) >= a->input->goodQuality) return false;

    Entity* edges[3];
    a->mesh->getDownward(worstTriangle, 1, edges);
    double longestLength = a->sizeField->measure(edges[0]);
    Entity* longestEdge = edges[0];
    for (int i=1; i<3; i++) {
      double length = a->sizeField->measure(edges[i]);
      if (length > longestLength) {
        longestLength = length;
        longestEdge = edges[i];
      }
    }

    Entity* oppositeVert = getTriVertOppositeEdge(a->mesh, worstTriangle, longestEdge);
    if (edgeSwap->run(longestEdge)) return true;
    if (splitCollapse.run(longestEdge, oppositeVert)) return true;
    for (int i=0; i<3; i++)
      if (edges[i] != longestEdge && collapseEdge(edges[i]))
        return true;
    return false;
  }

  bool isTwoLargeAngles(Entity* tet, Entity* problemEnts[4])
  {
    Entity* verts[4];
    a->mesh->getDownward(tet, 0, verts);
    Entity* vert = verts[0];
    Entity* face = getTetFaceOppositeVert(a->mesh, tet, vert);
    double area[4];
    int bit = getTetStats(a, vert, face, tet, problemEnts, area);
    return bit==3 || bit==5 || bit==6;
  }

  void fixTwoLargeAngles(Entity* tet, Entity* problemEnts[4])
  {
    if (edgeSwap->run(problemEnts[0])) return;
    if (edgeSwap->run(problemEnts[1])) return;
    Entity* v0 = getTriVertOppositeEdge(a->mesh, problemEnts[2], problemEnts[0]);
    if (splitCollapse.run(problemEnts[0], v0)) return;
    Entity* v1 = getTriVertOppositeEdge(a->mesh, problemEnts[3], problemEnts[1]);
    if (splitCollapse.run(problemEnts[1], v1)) return;
    if (doubleSplitCollapse.run(problemEnts)) return;
  }

  void fixThreeLargeAngles(Entity* tet, Entity* problemEnts[4])
  {
    Entity* faceEdges[3];
    a->mesh->getDownward(problemEnts[0], 1, faceEdges);
    if (edgeSwap->run(faceEdges[0])) return;
    if (edgeSwap->run(faceEdges[1])) return;
    if (edgeSwap->run(faceEdges[2])) return;
    if (runFaceSwap(a, problemEnts[0], true)) return;
    Entity* v = getTriVertOppositeEdge(a->mesh, problemEnts[2], problemEnts[1]);
    if (splitCollapse.run(problemEnts[1], v)) return;
    if (faceSplitCollapse.run(problemEnts[0], tet)) return;
  }

  void loopOnce()
  {
    Entity* tet;
    Iterator* it = a->mesh->begin(3);
    while ((tet = a->mesh->iterate(it))) {
      if (!getFlag(a, tet, BAD_QUALITY)) continue;
      clearFlag(a, tet, BAD_QUALITY);

      if (shortEdgeCase(tet)) continue;
      if (oneLargeAngle(tet)) continue;
      Entity* problemEnts[4];
      if (isTwoLargeAngles(tet, problemEnts))
        fixTwoLargeAngles(tet, problemEnts);
      else
        fixThreeLargeAngles(tet, problemEnts);
    }
  }
};

void fixElementShapesNew(Adapt* a)
{
  if ( ! a->input->shouldFixShape)
    return;
  double t0 = pcu::Time();
  int count = markBadQualityNew(a);
  print(a->mesh->getPCU(), "--iter %d of shape correction loop: #bad elements %d", 0, count);
  FixShape fixShape(a);
  int originalCount = count;
  int prev_count;
  int iter = 0;
  do {
    if ( ! count)
      break;
    prev_count = count;
    fixShape.loopOnce();
    if (a->mesh->getDimension() == 3)
      snap(a);
    count = markBadQualityNew(a);
    if (count >= prev_count)
      unMarkBadQualityNew(a); // to make sure markEntities does not complain!
    midBalance(a); // balance the mesh to avoid empty parts
    iter++;
    print(a->mesh->getPCU(), "loop %d: bad shapes went from %d to %d", iter, prev_count, count); 
  } while(count < prev_count);
  double t1 = pcu::Time();
  print(a->mesh->getPCU(), "bad shapes down from %d to %d in %f seconds", 
        originalCount,count,t1-t0);
}

}

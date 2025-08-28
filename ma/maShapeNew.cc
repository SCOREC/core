
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

static bool shortEdgeCase(Adapt* a, Entity* tet, Collapse& collapse)
{
  Entity* edges[6];
  a->mesh->getDownward(tet, 1, edges);
  bool shortEdgeCase = false;
  for (int i=0; i<6; i++) {
    if (a->sizeField->measure(edges[i]) < MINLENGTH) {
      shortEdgeCase = true;
      if (collapse.setEdge(edges[i]) && 
            collapse.checkClass() &&
            collapse.checkTopo() &&
            collapse.tryBothDirections(a->input->validQuality)) {
        collapse.destroyOldElements();
        return true;
      }
    }
  }
  return shortEdgeCase;
}

static double getWorstTriangle(Adapt* a, Entity* tet, Entity*& worstTriangle)
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

static bool oneLargeAngle(Adapt* a, Entity* tet, SingleSplitCollapse& splitCollapse, EdgeSwap* edgeSwap)
{
  Entity* worstTriangle;
  if (getWorstTriangle(a, tet, worstTriangle) >= a->input->goodQuality) return false;

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
  if (splitCollapse.run(longestEdge, oppositeVert))
    return true;
  if (edgeSwap->run(longestEdge))
    return true;
  return true;
}

static void fixLargeAngleTetsNew(Adapt* a)
{
  Collapse collapse;
  collapse.Init(a);
  FaceSplitCollapse faceSplitCollapse(a);
  SingleSplitCollapse splitCollapse(a);
  DoubleSplitCollapse doubleSplitCollapse(a);
  EdgeSwap* edgeSwap = makeEdgeSwap(a);

  Entity* tet;
  Iterator* it = a->mesh->begin(3);
  while ((tet = a->mesh->iterate(it))) {
    if (!getFlag(a, tet, BAD_QUALITY)) continue;
    clearFlag(a, tet, BAD_QUALITY);

    if (shortEdgeCase(a, tet, collapse)) continue;
    if (oneLargeAngle(a, tet, splitCollapse, edgeSwap)) continue;

    Entity* verts[4];
    a->mesh->getDownward(tet, 0, verts);
    Entity* vert = verts[0];
    Entity* face = getTetFaceOppositeVert(a->mesh, tet, vert);
    Entity* ents[4];
    double area[4];
    int bit = getTetStats(a, vert, face, tet, ents, area);

    if (bit==3 || bit==5 || bit==6) { //Two Large Angles
      if (doubleSplitCollapse.run(ents)) continue;
      if (edgeSwap->run(ents[0])) continue;
      if (edgeSwap->run(ents[1])) continue;
    }
    else { //Three Large Angles
      if (faceSplitCollapse.run(ents[0], tet)) continue;
      if (runFaceSwap(a, ents[0], true)) continue;
    }
  }
}

void fixElementShapesNew(Adapt* a)
{
  if ( ! a->input->shouldFixShape)
    return;
  double t0 = pcu::Time();
  int count = markBadQualityNew(a);
  print(a->mesh->getPCU(), "--iter %d of shape correction loop: #bad elements %d", 0, count);
  int originalCount = count;
  int prev_count;
  int iter = 0;
  do {
    if ( ! count)
      break;
    prev_count = count;
    fixLargeAngleTetsNew(a);
    /* We need to snap the new verts as soon as they are
     * created (to avoid future problems). At the moment
     * new verts are created only during 3D mesh adapt, so
     * we only run a bulk snap operation if the mesh is 3D.
     */
    if (a->mesh->getDimension() == 3)
      snap(a);
    count = markBadQualityNew(a);
    if (count >= prev_count)
      unMarkBadQualityNew(a); // to make sure markEntities does not complain!
    // balance the mesh to avoid empty parts
    midBalance(a);
    iter++;
  } while(count < prev_count);
  double t1 = pcu::Time();
  print(a->mesh->getPCU(), "bad shapes down from %d to %d in %f seconds", 
        originalCount,count,t1-t0);
}

}

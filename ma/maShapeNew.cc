
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

static void fixLargeAngleTetsNew(Adapt* a)
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

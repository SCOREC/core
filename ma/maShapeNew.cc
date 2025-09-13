
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

  int numCollapse=0;
  int numEdgeSwap=0;
  int numFaceSwap=0;
  int numEdgeSplitCollapse=0;
  int numFaceSplitCollapse=0;
  int numDoubleSplitCollapse=0;

  int numOneShortEdge=0;
  int numTwoShortEdge=0;
  int numThreeShortEdge=0;
  int numMoreShortEdge=0;
  int numOneLargeAngle=0;
  int numTwoLargeAngles=0;
  int numThreeLargeAngles=0;

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

  bool collapseToAdjacent(Entity* edge)
  {
    Entity* verts[2];
    a->mesh->getDownward(edge, 0, verts);
    for (int v=0; v<2; v++) {
      apf::Up adjEdges;
      a->mesh->getUp(verts[v], adjEdges);
      for (int e=0; e<adjEdges.n; e++) {
        if (adjEdges.e[e] == edge) continue;
        Entity* keep = getEdgeVertOppositeVert(a->mesh, adjEdges.e[e], verts[v]);
        bool alreadyFlagged = getFlag(a, keep, DONT_COLLAPSE);
        setFlag(a, keep, DONT_COLLAPSE);
        bool success = collapseEdge(adjEdges.e[e]);
        if (!alreadyFlagged) clearFlag(a, keep, DONT_COLLAPSE);
        if (success) return true;
      }
    }
    return false;
  }

  bool fixShortEdge(Entity* tet)
  {
    Entity* edges[6];
    a->mesh->getDownward(tet, 1, edges);
    for (int i=0; i<6; i++)
      if (a->sizeField->measure(edges[i]) < MINLENGTH)
        if (collapseEdge(edges[i]))
          { numCollapse++; return true; }
    for (int i=0; i<6; i++)
      if (a->sizeField->measure(edges[i]) < MINLENGTH)
        if (collapseToAdjacent(edges[i]))
          { numCollapse++; return true; }
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

  bool fixOneLargeAngle(Entity* tet)
  {
    Entity* worstTriangle;
    if (getWorstTriangle(tet, worstTriangle) >= .09) return false;

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
    if (edgeSwap->run(longestEdge))
      {numEdgeSwap++; return true;}
    if (splitCollapse.run(longestEdge, oppositeVert))
      {numEdgeSplitCollapse++; return true;}
    for (int i=0; i<3; i++)
      if (edges[i] != longestEdge && collapseEdge(edges[i]))
        {numCollapse++; return true;}
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
    if (edgeSwap->run(problemEnts[0])) {numEdgeSwap++; return;}
    if (edgeSwap->run(problemEnts[1])) {numEdgeSwap++; return;}
    Entity* faces[4];
    a->mesh->getDownward(tet, 2, faces);
    for (int i=0; i<4; i++) {
      if (isLowInHigh(a->mesh, faces[i], problemEnts[0])) {
        Entity* v0 = getTriVertOppositeEdge(a->mesh, faces[i], problemEnts[0]);
        if (splitCollapse.run(problemEnts[0], v0)) {numEdgeSplitCollapse++; return;}
      }
      else {
        Entity* v0 = getTriVertOppositeEdge(a->mesh, faces[i], problemEnts[1]);
        if (splitCollapse.run(problemEnts[1], v0)) {numEdgeSplitCollapse++; return;}
      }
    }
    if (doubleSplitCollapse.run(problemEnts)) {numDoubleSplitCollapse++; return;}
  }

  void fixThreeLargeAngles(Entity* tet, Entity* problemEnts[4])
  {
    Entity* tetEdges[6];
    a->mesh->getDownward(tet, 1, tetEdges);
    for (int i=0; i<6; i++)
      if (!isLowInHigh(a->mesh, problemEnts[0], tetEdges[i])) {
        Entity* verts[2];
        a->mesh->getDownward(tetEdges[i], 0, verts);
        Entity* keep = isLowInHigh(a->mesh, problemEnts[0], verts[0]) ? verts[0] : verts[1];
        bool alreadyFlagged = getFlag(a, keep, DONT_COLLAPSE);
        setFlag(a, keep, DONT_COLLAPSE);
        bool success = collapseEdge(tetEdges[i]);
        if (!alreadyFlagged) clearFlag(a, keep, DONT_COLLAPSE);
        if (success) {numCollapse++; return;}
      }
    Entity* faceEdges[3];
    a->mesh->getDownward(problemEnts[0], 1, faceEdges);
    if (edgeSwap->run(faceEdges[0])) {numEdgeSwap++; return;}
    if (edgeSwap->run(faceEdges[1])) {numEdgeSwap++; return;}
    if (edgeSwap->run(faceEdges[2])) {numEdgeSwap++; return;}
    if (runFaceSwap(a, problemEnts[0], true)) {numFaceSwap++; return;}
    Entity* v = getTriVertOppositeEdge(a->mesh, problemEnts[2], problemEnts[1]);
    if (splitCollapse.run(problemEnts[1], v)) {numEdgeSplitCollapse++; return;}
    if (faceSplitCollapse.run(problemEnts[0], tet)) {numFaceSplitCollapse++; return;}
  }

  void loopOnce()
  {
    Entity* tet;
    Iterator* it = a->mesh->begin(3);
    while ((tet = a->mesh->iterate(it))) {
      if (!getFlag(a, tet, BAD_QUALITY)) continue;
      clearFlag(a, tet, BAD_QUALITY);

      if (fixShortEdge(tet)) continue;
      if (fixOneLargeAngle(tet)) continue;
      Entity* problemEnts[4];
      if (isTwoLargeAngles(tet, problemEnts))
        fixTwoLargeAngles(tet, problemEnts);
      else
        fixThreeLargeAngles(tet, problemEnts);
    }
  }

  void printNumOperations()
  {
    print(a->mesh->getPCU(), "shape operations: \n collapses %17d\n edge swaps %16d\n face swaps %12d\n "
                              "edge split collapses %5d\n face split collapses %3d\n double split collapse %d\n",
                              numCollapse, numEdgeSwap, numFaceSwap, 
                              numEdgeSplitCollapse, numFaceSplitCollapse, numDoubleSplitCollapse);
  }

  bool isShortEdge(Entity* tet)
  {
    Entity* edges[6];
    a->mesh->getDownward(tet, 1, edges);
    int numShortEdges=0;
    for (int i=0; i<6; i++)
      if (a->sizeField->measure(edges[i]) < MINLENGTH)
        numShortEdges++;

    if (numShortEdges==0) return false;
    else if (numShortEdges==1) numOneShortEdge++;
    else if (numShortEdges==2) numTwoShortEdge++;
    else if (numShortEdges==3) numThreeShortEdge++;
    else numMoreShortEdge++;
    return true;
  }

  bool isOneLargeAngle(Entity* tet)
  {
    Entity* worstTriangle;
    if (getWorstTriangle(tet, worstTriangle) >= a->input->goodQuality) return false;
    return true;
  }

  void resetCounters()
  {
    numOneShortEdge=numTwoShortEdge=numThreeShortEdge=numMoreShortEdge=0;
    numOneLargeAngle=numTwoLargeAngles=numThreeLargeAngles=0;
  }

  void printBadShape(Entity* tet)
  {
    return;
    if (cbrt(a->shape->getQuality(tet)) > .11) return;
    EntityArray bad;
    bad.append(tet);
    apf::writeVtkFiles("shape_mesh", a->mesh);
    ma_dbg::createCavityMesh(a, bad, "shape_worst");

    Entity* faces[4];
    a->mesh->getDownward(tet, 2, faces);
    for (int f=0; f<4; f++) {
      apf::Up upward;
      a->mesh->getUp(faces[f], upward);
      for (int u=0; u<upward.n; u++)
        bad.append(upward.e[u]);
    }
    ma_dbg::createCavityMesh(a, bad, "shape_adjacent");

    Entity* problemEnts[4];
    if (isTwoLargeAngles(tet, problemEnts)) {
      printf("TWO LARGE ANGLES\n");
      setFlag(a, problemEnts[0], CHECKED);
      setFlag(a, problemEnts[1], CHECKED);
      ma_dbg::dumpMeshWithFlag(a, 0, 1, CHECKED, "shape_edges", "shape_edges");
    }
  }

  void printNumTypes()
  {
    resetCounters();
    Entity* tet;
    Iterator* it = a->mesh->begin(3);
    while ((tet = a->mesh->iterate(it))) {
      if (!getFlag(a, tet, BAD_QUALITY)) continue;
      printBadShape(tet);
      
      Entity* problemEnts[4];
      if (isShortEdge(tet)) continue;
      else if (isOneLargeAngle(tet))
        numOneLargeAngle++;
      else if (isTwoLargeAngles(tet, problemEnts))
        numTwoLargeAngles++;
      else
        numThreeLargeAngles++;
    }
    print(a->mesh->getPCU(), "bad shape types: \n oneShortEdge   \t%d\n twoShortEdges   \t%d\n threeShortEdges \t%d\n "
                              "moreShortEdges \t%d\n oneLargeAngle   \t%d\n twoLargeAngle   \t%d\n threeLargeAngle \t%d\n",
                              numOneShortEdge, numTwoShortEdge, numThreeShortEdge,
                              numMoreShortEdge, numOneLargeAngle, numTwoLargeAngles, numThreeLargeAngles);
  }
};

void fixElementShapesNew(Adapt* a)
{
  if ( ! a->input->shouldFixShape)
    return;
  double t0 = pcu::Time();
  bool oldForce = a->input->shouldForceAdaptation;
  a->input->shouldForceAdaptation = false;
  int count = markBadQualityNew(a);
  print(a->mesh->getPCU(), "loop %d: of shape correction loop: #bad elements %d", 0, count);
  FixShape fixShape(a);
  // fixShape.printNumTypes();
  int originalCount = count;
  int prev_count;
  int iter = 0;
  do {
    if (!count) break;
    prev_count = count;
    fixShape.loopOnce();
    if (a->mesh->getDimension() == 3)
      snap(a);
    count = markBadQualityNew(a);
    midBalance(a); // balance the mesh to avoid empty parts
    iter++;
    print(a->mesh->getPCU(), "loop %d: bad shapes went from %d to %d", iter, prev_count, count); 
  } while(count < prev_count);
  a->input->shouldForceAdaptation = oldForce;
  double t1 = pcu::Time();
  print(a->mesh->getPCU(), "bad shapes down from %d to %d in %f seconds", 
        originalCount,count,t1-t0);
  fixShape.printNumOperations();
  fixShape.printNumTypes();
  unMarkBadQualityNew(a);
}

}

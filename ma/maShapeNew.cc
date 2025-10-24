
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

class FixShape : public Operator
{
  public:
  Adapt* a;
  Mesh* mesh;
  Collapse collapse;
  SingleSplitCollapse splitCollapse;
  DoubleSplitCollapse doubleSplitCollapse;
  FaceSplitCollapse faceSplitCollapse;
  EdgeSwap* edgeSwap;
  Splits split;
  Entity* badTet;

  int numCollapse=0;
  int numEdgeSwap=0;
  int numFaceSwap=0;
  int numSplitReposition=0;
  int numRegionCollapse=0;
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

  FixShape(Adapt* adapt) : splitCollapse(adapt), doubleSplitCollapse(adapt), faceSplitCollapse(adapt), split(adapt)
  {
    a = adapt;
    mesh = a->mesh;
    collapse.Init(a);
    edgeSwap = makeEdgeSwap(a);
  }

  int getTargetDimension() { return 3; }
  bool shouldApply(Entity* e) {
    badTet = e;
    return getFlag(a, e, BAD_QUALITY); 
  }
  bool requestLocality(apf::CavityOp* o)
  {
    Entity* verts[4];
    a->mesh->getDownward(badTet, 0, verts);
    return o->requestLocality(verts, 4);
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
        if (!mesh->isOwned(keep)) continue; // TODO: improvement to requestLocality should remove use for this function
        if (mesh->isShared(keep)) continue; // TODO: improvement to requestLocality should remove use for this function
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

  double getWorstShape(EntityArray& tets, Entity*& worst)
  {
    double worstQuality = 1;
    for (int i=0; i<tets.size(); i++) {
      double quality = a->sizeField->measure(tets[i]);
      if (quality < worstQuality) {
        worstQuality = quality;
        worst = tets[i];
      }
    }
    return worstQuality;
  }

  Vector avgCavityPos(Entity* vert)
  {
    apf::Up edges;
    mesh->getUp(vert, edges);
    Vector avg(0,0,0);
    for (int i=0; i<edges.n; i++) {
      Entity* opp = getEdgeVertOppositeVert(mesh, edges.e[i], vert);
      avg += getPosition(mesh, opp);
    }
    avg = avg / edges.n;
    return avg;
  }

  void repositionVertex(Entity* vert)
  {
    EntityArray adjacent;
    mesh->getAdjacent(vert, mesh->getDimension(), adjacent);
    Entity* worstShape;
    Vector prevPos = getPosition(mesh, vert);
    getWorstShape(adjacent, worstShape);
    Vector dir = (avgCavityPos(vert) - prevPos).normalize();
    double speed = (avgCavityPos(vert) - prevPos).getLength()/2;
    double prevWorstQuality = 0;
    for (int i=0; i<10; i++) {
      double worstQuality = getWorstShape(adjacent, worstShape);
      if (worstQuality > prevWorstQuality) {
        prevWorstQuality = worstQuality;
        prevPos = getPosition(mesh, vert);
        mesh->setPoint(vert, 0, prevPos + (dir * speed));
        speed *= 2;
      }
      else {
        mesh->setPoint(vert, 0, prevPos);
        speed /= 4;
      }
    }
  }

  bool splitReposition(Entity* edge)
  {
    if (mesh->getModelType(mesh->toModel(edge)) != 3) return false;
    if (!split.setEdges(&edge, 1))
      return false;
    double worstQuality = getWorstQuality(a,split.getTets());
    split.makeNewElements();
    split.transfer();
    Entity* newVert = split.getSplitVert(0);
    repositionVertex(newVert);

    EntityArray adjacent;
    mesh->getAdjacent(newVert, mesh->getDimension(), adjacent);
    if (hasWorseQuality(a,adjacent,worstQuality)) {
      split.cancel();
      return false;
    }
    split.destroyOldElements();
    return true;
  }

  bool isOneLargeAngle(Entity* tet, Entity*& worstTriangle)
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
    return worstQuality < .09;
  }

  Entity* getLongestEdge(Entity* edges[3])
  {
    double longestLength = a->sizeField->measure(edges[0]);
    Entity* longestEdge = edges[0];
    for (int i=1; i<3; i++) {
      double length = a->sizeField->measure(edges[i]);
      if (length > longestLength) {
        longestLength = length;
        longestEdge = edges[i];
      }
    }
    return longestEdge;
  }

  bool fixOneLargeAngle(Entity* tet)
  {
    Entity* worstTriangle;
    if (!isOneLargeAngle(tet, worstTriangle)) return false; 
    Entity* edges[3];
    a->mesh->getDownward(worstTriangle, 1, edges);
    Entity* longestEdge = getLongestEdge(edges);
    Entity* oppositeVert = getTriVertOppositeEdge(a->mesh, worstTriangle, longestEdge);
    if (edgeSwap->run(longestEdge)) {numEdgeSwap++; return true;}
    if (splitCollapse.run(longestEdge, oppositeVert)) {numEdgeSplitCollapse++; return true;}
    for (int i=0; i<3; i++)
      if (edges[i] != longestEdge && collapseEdge(edges[i]))
        {numCollapse++; return true;}
    return false;
  }

  bool collapseRegion(Entity* tet, Entity* problemEnts[4]) 
  {
    Entity* faces[4];
    Entity* surface[4];
    Entity* interior[4];
    mesh->getDownward(tet, 2, faces);
    int s = 0, i = 0;
    for (int f=0; f<4; f++) {
      if (mesh->getModelType(mesh->toModel(faces[f])) == 2)
        surface[s++] = faces[f];
      else interior[i++] = faces[f];
    }
    if (s != 2 || i != 2) return false;
    Entity* in
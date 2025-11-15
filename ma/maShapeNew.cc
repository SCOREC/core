
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

FixShape::FixShape(Adapt* adapt) : splitCollapse(adapt), doubleSplitCollapse(adapt), faceSplitCollapse(adapt), split(adapt)
{
  a = adapt;
  mesh = a->mesh;
  collapse.Init(a);
  edgeSwap = makeEdgeSwap(a);
}

int FixShape::getTargetDimension() { return 3; }
bool FixShape::shouldApply(Entity* e) {
  badTet = e;
  return getFlag(a, e, BAD_QUALITY); 
}
bool FixShape::requestLocality(apf::CavityOp* o)
{
  Entity* verts[4];
  a->mesh->getDownward(badTet, 0, verts);
  return o->requestLocality(verts, 4);
}

bool FixShape::collapseEdge(Entity* edge)
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

bool FixShape::collapseToAdjacent(Entity* edge)
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

bool FixShape::fixShortEdge(Entity* tet)
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

double FixShape::getWorstShape(EntityArray& tets, Entity*& worst)
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

Vector FixShape::avgCavityPos(Entity* vert)
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

void FixShape::repositionVertex(Entity* vert)
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

bool FixShape::splitReposition(Entity* edge)
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

bool FixShape::isOneLargeAngle(Entity* tet, Entity*& worstTriangle)
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

Entity* FixShape::getLongestEdge(Entity* edges[3])
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

bool FixShape::fixOneLargeAngle(Entity* tet)
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

bool FixShape::collapseRegion(Entity* tet, Entity* problemEnts[4]) 
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
  Entity* interiorEdge = mesh->getModelType(mesh->toModel(problemEnts[0])) == 3 ? problemEnts[0] : problemEnts[1];
  Entity* surfaceEdge = problemEnts[0] == interiorEdge ? problemEnts[1] : problemEnts[0];
  if (mesh->getModelType(mesh->toModel(surfaceEdge)) == 1)
    return false;
  if (!isLowInHigh(mesh, surface[0], surfaceEdge) || !isLowInHigh(mesh, surface[1], surfaceEdge))
    return false;
  Model* modelFace = mesh->toModel(surface[0]);
  mesh->destroy(tet);
  mesh->destroy(surface[0]);
  mesh->destroy(surface[1]);
  mesh->destroy(surfaceEdge);
  
  mesh->setModelEntity(interior[0], modelFace);
  mesh->setModelEntity(interior[1], modelFace);
  mesh->setModelEntity(interiorEdge, modelFace);
  return true;
}

bool FixShape::isTwoLargeAngles(Entity* tet, Entity* problemEnts[4])
{
  Entity* verts[4];
  a->mesh->getDownward(tet, 0, verts);
  Entity* vert = verts[0];
  Entity* face = getTetFaceOppositeVert(a->mesh, tet, vert);
  double area[4];
  int bit = getTetStats(a, vert, face, tet, problemEnts, area);
  return bit==3 || bit==5 || bit==6;
}

bool FixShape::fixTwoLargeAngles(Entity* tet, Entity* problemEnts[4])
{
  if (collapseRegion(tet, problemEnts)) {numRegionCollapse++; return true;}
  if (edgeSwap->run(problemEnts[0])) {numEdgeSwap++; return true;}
  if (edgeSwap->run(problemEnts[1])) {numEdgeSwap++; return true;}
  if (splitReposition(problemEnts[0])) {numSplitReposition++; return true;}
  if (splitReposition(problemEnts[1])) {numSplitReposition++; return true;}
  Entity* faces[4];
  a->mesh->getDownward(tet, 2, faces);
  for (int i=0; i<4; i++) {
    if (isLowInHigh(a->mesh, faces[i], problemEnts[0])) {
      Entity* v0 = getTriVertOppositeEdge(a->mesh, faces[i], problemEnts[0]);
      if (splitCollapse.run(problemEnts[0], v0)) {numEdgeSplitCollapse++; return true;}
    }
    else {
      Entity* v0 = getTriVertOppositeEdge(a->mesh, faces[i], problemEnts[1]);
      if (splitCollapse.run(problemEnts[1], v0)) {numEdgeSplitCollapse++; return true;}
    }
  }
  if (doubleSplitCollapse.run(problemEnts)) {numDoubleSplitCollapse++; return true;}
  return false;
}

bool FixShape::fixThreeLargeAngles(Entity* tet, Entity* problemEnts[4])
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
      if (success) {numCollapse++; return true;}
    }
  Entity* faceEdges[3];
  a->mesh->getDownward(problemEnts[0], 1, faceEdges);
  if (edgeSwap->run(faceEdges[0])) {numEdgeSwap++; return true;}
  if (edgeSwap->run(faceEdges[1])) {numEdgeSwap++; return true;}
  if (edgeSwap->run(faceEdges[2])) {numEdgeSwap++; return true;}
  if (runFaceSwap(a, problemEnts[0], true)) {numFaceSwap++; return true;}
  Entity* v = getTriVertOppositeEdge(a->mesh, problemEnts[2], problemEnts[1]);
  if (splitCollapse.run(problemEnts[1], v)) {numEdgeSplitCollapse++; return true;}
  if (faceSplitCollapse.run(problemEnts[0], tet)) {numFaceSplitCollapse++; return true;}
  return false;
}

void FixShape::apply()
{
  clearFlag(a, badTet, BAD_QUALITY);
  if (fixShortEdge(badTet)) return;
  if (fixOneLargeAngle(badTet)) return;
  Entity* problemEnts[4];
  if (isTwoLargeAngles(badTet, problemEnts))
    fixTwoLargeAngles(badTet, problemEnts);
  else
    fixThreeLargeAngles(badTet, problemEnts);
}

void FixShape::resetCounters()
{
  numOneShortEdge=numTwoShortEdge=numThreeShortEdge=numMoreShortEdge=0;
  numOneLargeAngle=numTwoLargeAngles=numThreeLargeAngles=0;
}
int FixShape::collect(int val) {
  return a->mesh->getPCU()->Add<long>(val);
}
void FixShape::printNumOperations()
{
  print(a->mesh->getPCU(), "shape operations: \n collapses %17d\n edge swaps %16d\n split reposition %9d\n double split collapse %d\n "
                            "edge split collapses %5d\n face split collapses %3d\n region collapses %7d\n face swaps %12d\n ",
                            collect(numCollapse), collect(numEdgeSwap), collect(numSplitReposition), collect(numDoubleSplitCollapse),
                            collect(numEdgeSplitCollapse), collect(numFaceSplitCollapse), collect(numRegionCollapse), collect(numFaceSwap));
}
void FixShape::printBadTypes()
{
  print(a->mesh->getPCU(), "bad shape types: \n oneShortEdge   \t%d\n twoShortEdges   \t%d\n threeShortEdges \t%d\n "
                            "moreShortEdges \t%d\n oneLargeAngle   \t%d\n twoLargeAngle   \t%d\n threeLargeAngle \t%d\n",
                            collect(numOneShortEdge), collect(numTwoShortEdge), collect(numThreeShortEdge),
                            collect(numMoreShortEdge), collect(numOneLargeAngle), collect(numTwoLargeAngles), collect(numThreeLargeAngles));
}

bool FixShape::isShortEdge(Entity* tet)
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

void FixShape::printBadShape(Entity* problemTet)
{
  Entity* worstTri;
  Entity* problemEnts[4];
  if (isShortEdge(problemTet)) print(a->mesh->getPCU(), "Worst is short\n");
  else if (isOneLargeAngle(problemTet, worstTri)) print(a->mesh->getPCU(), "Worst is one large angle\n");
  else if (isTwoLargeAngles(problemTet, problemEnts)) print(a->mesh->getPCU(), "Worst is two large angles\n");
  else print(a->mesh->getPCU(), "Worst is three large angles\n");

  ma_dbg::addClassification(a);

  EntitySet bad;
  bad.insert(problemTet);
  EntitySet adjacent1 = getNextLayer(a, bad);
  EntitySet adjacent2 = getNextLayer(a, adjacent1);

  ma_dbg::flagEntityAllDim(a, 3, "worst_tet", &problemTet, 1);
  ma_dbg::flagEntity(a, 3, "worst_tet", bad);
  ma_dbg::flagEntity(a, 3, "worst_tet_adj", adjacent1);
  ma_dbg::flagEntity(a, 3, "worst_tet_adj2", adjacent2);

  std::vector<Entity*> badFaces;
  Iterator* it = a->mesh->begin(3);
  Entity* tet;
  while ((tet = a->mesh->iterate(it))) {
    if (!getFlag(a, tet, BAD_QUALITY)) continue;
    Entity* faces[4];
    mesh->getDownward(tet, 2, faces);
    for (Entity* face : faces)
      if (mesh->getModelType(mesh->toModel(face)) == 2)
        badFaces.push_back(face);
  }
  ma_dbg::flagEntity(a, 2, "bad_surface_tets", &badFaces[0], badFaces.size());
  apf::writeVtkFiles("mesh_tets", a->mesh, 3);
  apf::writeVtkFiles("mesh_faces", a->mesh, 2);
  apf::writeVtkFiles("mesh_edges", a->mesh, 1);
}

void FixShape::printNumTypes()
{
  resetCounters();
  Entity* tet;
  Iterator* it = a->mesh->begin(3);
  double worstQual = 1;
  Entity* worstShape;
  while ((tet = a->mesh->iterate(it))) {
    if (!getFlag(a, tet, BAD_QUALITY)) continue;
    double qual = a->shape->getQuality(tet);
    if (qual < worstQual) {
      worstQual = qual;
      worstShape = tet;
    }
    
    Entity* problemEnts[4];
    Entity* worst;
    if (isShortEdge(tet)) continue;
    else if (isOneLargeAngle(tet, worst))
      numOneLargeAngle++;
    else if (isTwoLargeAngles(tet, problemEnts))
      numTwoLargeAngles++;
    else
      numThreeLargeAngles++;
  }
  printBadTypes();
  // printBadShape(worstShape);
}

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
    applyOperator(a,&fixShape);
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
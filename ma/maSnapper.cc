/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maSnapper.h"
#include "maAdapt.h"
#include "maShapeHandler.h"
#include "maSnap.h"
#include "maDBG.h"
#include <apfCavityOp.h>
#include <pcu_util.h>
#include <lionPrint.h>
#include <iostream>

namespace ma {

Snapper::Snapper(Adapt* a, Tag* st, bool is)
{
  adapter = a;
  snapTag = st;
  collapse.Init(a);
  isSimple = is;
  dug = false;
  vert = 0;
}

void Snapper::setVert(Entity* v)
{
  vert = v;
}

Entity* Snapper::getVert()
{
  return vert;
}

bool Snapper::requestLocality(apf::CavityOp* o)
{
  if (!o->requestLocality(&vert, 1))
    return false;
/* in order to try an edge collapse (we don't yet know
   which edge), bring in a cavity such that all adjacent
   edges have both vertices local.
   This is basically two layers of elements around the vertex */
  apf::Up edges;
  adapter->mesh->getUp(vert,edges);
  apf::Up ovs;
  ovs.n = edges.n;
  for (int i = 0; i < edges.n; ++i)
    ovs.e[i] = apf::getEdgeVertOppositeVert(adapter->mesh, edges.e[i], vert);
  return o->requestLocality(&ovs.e[0], ovs.n);
}

static void computeNormals(Mesh* m, Upward& es, apf::NewArray<Vector>& normals)
{
  if (m->getDimension() != 2)
    return;
  normals.allocate(es.getSize());
  for (size_t i = 0; i < es.getSize(); ++i)
    normals[i] = getTriNormal(m, es[i]);
}

static bool didInvert(Mesh* m, Vector& oldNormal, Entity* tri)
{
  return (oldNormal * getTriNormal(m, tri)) < 0;
}

static void updateVertexParametricCoords(
    Mesh* m,
    Entity* vert,
    Vector& newTarget)
{
  PCU_ALWAYS_ASSERT_VERBOSE(m->getType(vert) == apf::Mesh::VERTEX,
      "expecting a vertex!");

  // if vert is classified on a model vert or edge return
  Model* g = m->toModel(vert);
  if (m->getModelType(g) != 2)
    return;

  // get the list of upward adj edges that are
  // classified on the same model face as vert
  apf::Up edges;
  m->getUp(vert,edges);
  apf::Up oes;
  oes.n = edges.n;
  int counter = 0;
  for (int i = 0; i < edges.n; ++i) {
    Model* h = m->toModel(edges.e[i]);
    if (m->getModelType(h) == 3)
      continue;
    PCU_ALWAYS_ASSERT_VERBOSE(g == h,
    	"expecting the model to be the same for current edge and vert");
    oes.e[counter] = edges.e[i];
    counter++;
  }

  Vector pBar(0., 0., 0.);
  for (int i = 0; i < counter; i++) {
    Vector pTmp;
    transferParametricOnEdgeSplit(m, oes.e[i], 0.5, pTmp);
    pBar += pTmp;
  }
  pBar = pBar / oes.n;

  m->snapToModel(m->toModel(vert), pBar, newTarget);
  m->setParam(vert, pBar);
}

static void printFPP(Adapt* a, FirstProblemPlane* FPP)
{
  apf::writeVtkFiles("FPP_Mesh", a->mesh);
  EntityArray invalid;
  for (int i=0; i<FPP->problemRegions.n; i++){
    invalid.append(FPP->problemRegions.e[i]);
  }
  ma_dbg::createCavityMesh(a, invalid, "FPP_Invalid");

  for (int i=0; i<FPP->commEdges.n; i++) setFlag(a, FPP->commEdges.e[i], CHECKED);
  ma_dbg::dumpMeshWithFlag(a, 0, 1, CHECKED, "FPP_CommEdges", "FPP_CommEdges");
  for (int i=0; i<FPP->commEdges.n; i++) clearFlag(a, FPP->commEdges.e[i], CHECKED);

  setFlag(a, FPP->vert, CHECKED);
  ma_dbg::dumpMeshWithFlag(a, 0, 0, CHECKED, "FPP_Vertex", "FPP_Vertex");
  clearFlag(a, FPP->vert, CHECKED);

  setFlag(a, FPP->problemFace, CHECKED);
  ma_dbg::dumpMeshWithFlag(a, 0, 2, CHECKED, "FPP_Face", "FPP_Face");
  clearFlag(a, FPP->problemFace, CHECKED);

  setFlag(a, FPP->problemRegion, CHECKED);
  ma_dbg::dumpMeshWithFlag(a, 0, 3, CHECKED, "FPP_Region", "FPP_Region");
  clearFlag(a, FPP->problemRegion, CHECKED);
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
      collapse.tryBothDirections(0)) {
    collapse.destroyOldElements();
    result = true;
  }  
  if (!alreadyFlagged) clearFlag(a, keep, DONT_COLLAPSE);
  return result;
}

struct BestCollapse
{
  double quality=-1;
  Entity* edge;
  Entity* keep;
};

static void getBestQualityCollapse(Adapt* a, Entity* edge, Entity* keep, Collapse& collapse, BestCollapse& best)
{
  PCU_ALWAYS_ASSERT(a->mesh->getType(edge) == apf::Mesh::EDGE);
  bool alreadyFlagged = true;
  if (keep) alreadyFlagged = getFlag(a, keep, DONT_COLLAPSE);
  if (!alreadyFlagged) setFlag(a, keep, DONT_COLLAPSE);
  if (collapse.setEdge(edge) &&
      collapse.checkClass() &&
      collapse.checkTopo()) {
    double quality = collapse.getQualityFromCollapse();
    if (quality > best.quality) {
      best.quality = quality;
      best.edge = edge;
      best.keep = keep;
    }
  }
  if (!alreadyFlagged) clearFlag(a, keep, DONT_COLLAPSE);
}

static bool sameSide(Adapt* a, Entity* testVert, Entity* refVert, Entity* face)
{
  Entity* faceVert[3];
  a->mesh->getDownward(face, 0, faceVert);
  Vector facePos[3];
  for (int i=0; i < 3; ++i)
    facePos[i] = getPosition(a->mesh,faceVert[i]);
  
  Vector normal = apf::cross((facePos[1]-facePos[0]),(facePos[2]-facePos[0]));
  Vector testPos = getPosition(a->mesh, testVert);
  Vector refPos = getPosition(a->mesh, refVert);
  const double tol=1e-12;

  double dr = (testPos - facePos[0]) * normal;
  if (dr*dr < tol) return false; //testVert is on the face
  double ds = (refPos - facePos[0]) * normal;
  if (dr*ds < 0.0) return false; //different sides of face
  return true; //same side of face
}

static bool tryCollapseTetEdges(Adapt* a, Collapse& collapse, FirstProblemPlane* FPP)
{
  apf::Up& commEdges = FPP->commEdges;
  BestCollapse best;

  for (int i=0; i<commEdges.n; i++) {
    getBestQualityCollapse(a, commEdges.e[i], 0, collapse, best);
  }

  for (int i=0; i<commEdges.n; i++) {
    Entity* edge = commEdges.e[i];
    Entity* vertexFPP = getEdgeVertOppositeVert(a->mesh, edge, FPP->vert);
    apf::Up adjEdges;
    a->mesh->getUp(vertexFPP, adjEdges);
    for (int j=0; j<adjEdges.n; j++) {
      Entity* edgeDel = adjEdges.e[j];
      if (edgeDel==edge) continue;
      if (isLowInHigh(a->mesh, FPP->problemFace, edgeDel)) continue;
      Entity* vertKeep = getEdgeVertOppositeVert(a->mesh, edgeDel, vertexFPP);
      if (sameSide(a, vertKeep, FPP->vert, FPP->problemFace)) continue;
      getBestQualityCollapse(a, edgeDel, vertKeep, collapse, best);
    }
  }

  if (best.quality > 0) 
    return tryCollapseEdge(a, best.edge, best.keep, collapse);
  else return false;
}

static bool tryReduceCommonEdges(Adapt* a, Collapse& collapse, FirstProblemPlane* FPP)
{
  apf::Up& commEdges = FPP->commEdges;
  BestCollapse best;

  Entity* pbEdges[3];
  a->mesh->getDownward(FPP->problemFace, 1, pbEdges);
  switch(commEdges.n) {
    case 2: {
      Entity* v1 = getEdgeVertOppositeVert(a->mesh, commEdges.e[0], FPP->vert);
      Entity* v2 = getEdgeVertOppositeVert(a->mesh, commEdges.e[1], FPP->vert);
      
      for (int i=0; i<3; i++) {
        Entity* pbVert[2];
        a->mesh->getDownward(pbEdges[i], 0, pbVert);
        if (pbVert[0] == v1 && pbVert[1] == v2) break;
        if (pbVert[1] == v1 && pbVert[0] == v2) break;
        getBestQualityCollapse(a, pbEdges[i], 0, collapse, best);
      }
      break;
    }
    case 3: {
      for (int i=0; i<3; i++)
        getBestQualityCollapse(a, pbEdges[i], 0, collapse, best);
      break;
    }
  }
  if (best.quality > 0) 
    return tryCollapseEdge(a, best.edge, best.keep, collapse);
  else return false;
}

static bool tryCollapseToVertex(Adapt* a, Collapse& collapse, FirstProblemPlane* FPP)
{
  Vector position = getPosition(a->mesh, FPP->vert);
  Vector target;
  a->mesh->getDoubleTag(FPP->vert, FPP->snapTag, &target[0]);
  double distTarget = (position - target).getLength();

  BestCollapse best;

  for (size_t i = 0; i < FPP->commEdges.n; ++i) {
    Entity* edge = FPP->commEdges.e[i];
    Entity* vertexOnFPP = getEdgeVertOppositeVert(a->mesh, edge, FPP->vert);
    Vector vFPPCoord = getPosition(a->mesh, vertexOnFPP);
    //TODO: add logic for boundary layers
    double distToFPPVert = (vFPPCoord - target).getLength();
    if (distToFPPVert > distTarget) continue;
    getBestQualityCollapse(a, edge, FPP->vert, collapse, best);
  }

  if (best.quality > 0) 
    return tryCollapseEdge(a, best.edge, best.keep, collapse);
  else return false;
}

static FirstProblemPlane* getFPP(Adapt* a, Entity* vertex, Tag* snapTag, apf::Up& invalid)
{
  FirstProblemPlane* FPP = new FirstProblemPlane(a, snapTag);
  FPP->setVertex(vertex);
  FPP->setBadElements(invalid);
  std::vector<Entity*> commEdges;
  FPP->getCandidateEdges(commEdges);
  return FPP;
}

static void getInvalidTets(Mesh* mesh, Upward& adjacentElements, apf::Up& invalid)
{
  invalid.n = 0;
  Vector v[4];
  for (size_t i = 0; i < adjacentElements.getSize(); ++i) {
    ma::getVertPoints(mesh,adjacentElements[i],v);
    if ((cross((v[1] - v[0]), (v[2] - v[0])) * (v[3] - v[0])) < 0)
      invalid.e[invalid.n++] = adjacentElements[i];
  }
}

static bool tryReposition(Adapt* adapt, Entity* vertex, Tag* snapTag, apf::Up& invalid) 
{
  Mesh* mesh = adapt->mesh;
  if (!mesh->hasTag(vertex, snapTag)) return true;
  Vector prev = getPosition(mesh, vertex);
  Vector target;
  mesh->getDoubleTag(vertex, snapTag, &target[0]);
  Upward adjacentElements;
  mesh->getAdjacent(vertex, mesh->getDimension(), adjacentElements);
  mesh->setPoint(vertex, 0, target);
  getInvalidTets(mesh, adjacentElements, invalid);
  if (invalid.n == 0) return true;
  mesh->setPoint(vertex, 0, prev);
  return false;
}

bool Snapper::trySimpleSnap()
{
  apf::Up invalid;
  return tryReposition(adapter, vert, snapTag, invalid);
}

static int numFailed = 0;

bool Snapper::run()
{
  apf::Up invalid;
  bool success = tryReposition(adapter, vert, snapTag, invalid);

  if (success) {
    adapter->mesh->removeTag(vert,snapTag);
    clearFlag(adapter, vert, SNAP);
    return true;
  }

  FirstProblemPlane* FPP = getFPP(adapter, vert, snapTag, invalid);

  if (!success) success = tryCollapseToVertex(adapter, collapse, FPP);
  // if (!success) success = tryReduceCommonEdges(a, collapse, FPP);
  // if (!success) success = tryCollapseTetEdges(adapter, collapse, FPP);

  if (!success) numFailed++;
  if (!success && numFailed == 1) printFPP(adapter, FPP);
  
  if (FPP) delete FPP;
  return !success;
}

FirstProblemPlane::FirstProblemPlane(Adapt* a, Tag* st)
{
  adapter = a;
  snapTag = st;
  problemFace = 0;
  problemRegion = 0;
  commEdges.n = 0;
  tol = 1.0e-14;
}

void FirstProblemPlane::setVertex(Entity* v)
{
  vert = v;
}

void FirstProblemPlane::setBadElements(apf::Up& badElements)
{
  problemRegions.n = badElements.n;
  for (int i = 0; i < badElements.n; i++) {
    problemRegions.e[i] = badElements.e[i];
  }
}


void FirstProblemPlane::getCandidateEdges(std::vector<Entity*> &edges)
{
  edges.clear();
  if (find())
    findCandidateEdges(edges);
}

bool FirstProblemPlane::find()
{
  Mesh* mesh = adapter->mesh;
  std::vector<double> dists;
  double minDist = 1.0e6;

  // determine distances to all possible problem faces, the shortest
  // distance and its intersection on first problem plane (FPP)
  int n;
  n = problemRegions.n;
  Entity* elem;
  Entity* face;
  Ray ray;
  Vector target;
  mesh->getDoubleTag(vert, snapTag, &target[0]);

  ray.start = getPosition(mesh, vert);
  ray.dir   = target - ray.start;

  dists.clear();
  for (int i = 0; i < n; i++) {
    elem = problemRegions.e[i];
    face = getTetFaceOppositeVert(mesh, elem, vert);
    std::vector<Vector> coords;
    getFaceCoords(mesh, face, coords);

    Vector intersect;
    bool isInf;
    bool ok = intersectRayFace(ray, coords, intersect, isInf);

    if (ok){
      if (isInf)
        lion_oprint(1, "Info: Found Infinitely Many Intersection Points!\n");
      Vector newDirection = intersect - ray.start;
      if (newDirection.getLength() < minDist) {
	      dists.push_back(newDirection.getLength());
      	minDist = dists.back();
      	problemFace = face;
      	problemRegion = elem;
      	intersection = intersect;
      	// do not need to check whether the move is valid since the valid
      	// ones should have been taken care of by this point
      }
      else
      	dists.push_back(newDirection.getLength());
    }
  }


  apf::Up coplanarProblemRegions;
  coplanarProblemRegions.n = 0;

  if (!problemRegion) {
    problemRegion = problemRegions.e[0];
    problemFace = getTetFaceOppositeVert(mesh, problemRegion, vert);
    coplanarProblemRegions.n = n;
    for (int i = 0; i < n; i++) {
      coplanarProblemRegions.e[i] = problemRegions.e[i];
    }
  }
  else {
    minDist += tol;
    for (int i = 0; i < n; i++) {
      if (dists[i] < minDist) {
        coplanarProblemRegions.e[coplanarProblemRegions.n] = problemRegions.e[i];
        coplanarProblemRegions.n++;
      }
    }
  }


  findCommonEdges(coplanarProblemRegions);
  return true;
}

void FirstProblemPlane::findCandidateEdges(std::vector<Entity*> &edges)
{
  edges.clear();
  Mesh* mesh = adapter->mesh;

  // We deny collapsing that moves further away form the current target
  // Need the original dist b/w the current vert and the target snap point
  Vector x, t;
  x = getPosition(mesh, vert);
  mesh->getDoubleTag(vert, snapTag, &t[0]);

  double dist = (x - t).getLength();

  Entity* edge;
  Entity* v;

  for (int i = 0; i < commEdges.n; i++) {
    edge = commEdges.e[i];
    Downward dv;
    mesh->getDownward(edge, 0, dv);
    (dv[0] == vert) ? v = dv[1] : v = dv[0];
    Vector vCoord = getPosition(mesh, v);

    if (false) {;} // boundary layer stuff

    double candidateDist = (vCoord - t).getLength();
    if (candidateDist > dist)
      continue;
    else
      edges.push_back(edge);
  }
}



bool
FirstProblemPlane::intersectRayFace(const Ray& ray, const std::vector<Vector>& coords,
    Vector& intersection, bool& isInf)
{
  bool res = false;
  isInf = false;
  if (coords.size() != 3){
    lion_oprint(1,"coords.size() is %d\n", coords.size());
    lion_oprint(1,"No implementation for non-tri faces!\n");
    res = false;
  }

  PCU_ALWAYS_ASSERT(ray.dir.getLength() > tol);
  Vector start = ray.start;
  Vector dir   = ray.dir;

  Vector p0p1 = coords[1] - coords[0];
  Vector p0p2 = coords[2] - coords[0];
  Vector startP0 = coords[0] - start;

  Vector faceAreaVect = apf::cross(p0p1, p0p2);
  double faceAreaSize = faceAreaVect.getLength();

  double vol = std::fabs(dir * faceAreaVect);
  double volPrime = std::fabs(startP0 * faceAreaVect);

  if (vol <= tol * faceAreaSize) { // dir _|_ face consisting of coords
    if (volPrime <= tol * faceAreaSize) {
      isInf = true;
      res = true;
      intersection = (coords[0] + coords[1] + coords[2]) * (1./3.);
    }
    else {
      res = false;
    }
  }
  else {
    intersection = start + dir * (volPrime / vol);
    Vector newDir = intersection - start;
    if (newDir * dir < 0)
    {
      res = false;
    }
    else
      res = true;
  }
  return res;
}

void FirstProblemPlane::findCommonEdges(apf::Up& cpRegions)
{
  Mesh* mesh = adapter->mesh;
  if (cpRegions.n == 1) {
    Downward edges;
    int nDownEdges = mesh->getDownward(cpRegions.e[0], 1, edges);
    for (int i = 0; i < nDownEdges; i++) {
      if (isLowInHigh(mesh, edges[i], vert)) {
        commEdges.e[commEdges.n] = edges[i];
        commEdges.n++;
      }
    }
    return;
  }
  // determine the problem face closest to intersection
  Entity* region;
  Entity* tmpRegion;
  Entity* face;

  Vector ctrToIntersect = getCenter(mesh, problemFace) - intersection;
  double minDist = ctrToIntersect.getLength();
  tmpRegion = problemRegion;

  for (int i = 0; i < cpRegions.n; i++) {
    region = cpRegions.e[i];
    if (region == tmpRegion) continue;
    face = getTetFaceOppositeVert(mesh, region, vert);
    ctrToIntersect = getCenter(mesh, face) - intersection;
    double dist = ctrToIntersect.getLength();
    if (dist < minDist) {
      minDist = dist;
      problemFace = face;
      problemRegion = region;
    }
  }

  Downward edges;
  int flag = 0;
  int nDownEdges = mesh->getDownward(problemRegion, 1, edges);
  for (int i = 0; i < nDownEdges; i++) {
    if (isLowInHigh(mesh, edges[i], vert)) {
      flag = 1;
      for (int j = 0; j < cpRegions.n; j++) {
      	region = cpRegions.e[j];
      	if (region == problemRegion) continue;
      	if (!isLowInHigh(mesh, region, edges[i])) {
      	  flag = 0;
      	  break;
        }
      }
      if (flag) {
        commEdges.e[commEdges.n] = edges[i];
        commEdges.n++;
      }
    }
  }
}

Entity* getTetFaceOppositeVert(Mesh* m, Entity* e, Entity* v)
{
  Downward faces;
  Entity* oppositeFace = 0;
  int nDownFaces = m->getDownward(e, 2, faces);
  for (int i = 0; i < nDownFaces; i++) {
    Downward verts;
    int nDownVerts = m->getDownward(faces[i], 0, verts);
    int j;
    for (j = 0; j < nDownVerts; j++) {
      if (v == verts[j])
      	break;
    }
    if (j == nDownVerts)
      oppositeFace = faces[i];
    else
      continue;
  }

  // make sure that oppositeFace is what it is meant to be!
  Downward verts;
  int numDownVerts = m->getDownward(oppositeFace, 0, verts);
  bool flag = true;
  for (int i = 0; i < numDownVerts; i++) {
    if (v == verts[i]){
      flag = false;
      break;
    }
  }

  PCU_ALWAYS_ASSERT(flag);

  return oppositeFace;
}

void getFaceCoords(Mesh* m, Entity* face, std::vector<Vector>& coords)
{
  Downward verts;
  int nDownVerts = m->getDownward(face, 0, verts);
  PCU_ALWAYS_ASSERT(nDownVerts);
  for (int i = 0; i < nDownVerts; i++)
    coords.push_back(getPosition(m, verts[i]));
}

Vector getCenter(Mesh* mesh, Entity* face)
{
  PCU_ALWAYS_ASSERT(face);
  Downward verts;
  int nDownVerts = mesh->getDownward(face, 0, verts);
  PCU_ALWAYS_ASSERT(nDownVerts == 3);
  Vector center(0., 0., 0.);
  for (int i = 0; i < nDownVerts; i++)
    center += getPosition(mesh, verts[i]);

  center = center / 3.;
  return center;
}

bool isLowInHigh(Mesh* mesh, Entity* highEnt, Entity* lowEnt)
{
  PCU_ALWAYS_ASSERT(mesh->getType(highEnt) > mesh->getType(lowEnt));
  Downward down;
  int nDown = mesh->getDownward(highEnt, apf::getDimension(mesh, lowEnt), down);
  for (int i = 0; i < nDown; i++) {
    if (lowEnt == down[i])
      return true;
  }
  return false;
}

}

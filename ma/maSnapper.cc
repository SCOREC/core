/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maSnapper.h"
#include "maAdapt.h"
#include "maShapeHandler.h"
#include <apfCavityOp.h>
#include <cassert>
#include <iostream>

namespace ma {

Snapper::Snapper(Adapt* a, Tag* st, bool is)
{
  adapter = a;
  snapTag = st;
  collapse.Init(a);
  isSimple = is;
  dug = false;
  toFPP = false;
  vert = 0;
}

bool Snapper::setVert(Entity* v, apf::CavityOp* o)
{
  vert = v;
  if (!o->requestLocality(&vert, 1))
    return false;
  if (isSimple)
    return true;
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

static void collectBadElements(Adapt* a, Upward& es,
    apf::NewArray<Vector>& normals, apf::Up& bad)
{
  Mesh* m = a->mesh;
  bad.n = 0;
  for (size_t i = 0; i < es.getSize(); ++i) {
/* for now, when snapping a vertex on the boundary
   layer, ignore the quality of layer elements.
   not only do we not have metrics for this, but the
   algorithm that moves curves would need to change */
    if (getFlag(a, es[i], LAYER))
      continue;
    double quality = a->shape->getQuality(es[i]);
    if (quality < a->input->validQuality)
      bad.e[bad.n++] = es[i];
/* check for triangles whose normals have changed by
   more than 90 degrees when the vertex was snapped */
    else if ((m->getDimension() == 2) &&
             didInvert(m, normals[i], es[i]))
      bad.e[bad.n++] = es[i];
  }
  assert(bad.n < (int)(sizeof(bad.e) / sizeof(Entity*)));
}

static bool trySnapping(Adapt* adapter, Tag* tag, Entity* vert,
    apf::Up& badElements)
{
  Mesh* mesh = adapter->mesh;
  Vector x = getPosition(mesh, vert);
  Vector s;
  mesh->getDoubleTag(vert, tag, &s[0]);
/* gather the adjacent elements */
  Upward elements;
  mesh->getAdjacent(vert, mesh->getDimension(), elements);
/* in 2D, get the old triangle normals */
  apf::NewArray<Vector> normals;
  computeNormals(mesh, elements, normals);
/* move the vertex to the desired point */
  mesh->setPoint(vert, 0, s);
/* check resulting cavity */
  collectBadElements(adapter, elements, normals, badElements);
  if (badElements.n) {
    /* not ok, put the vertex back where it was */
    mesh->setPoint(vert, 0, x);
    return false;
  } else {
    /* ok, take off the snap tag */
    mesh->removeTag(vert, tag);
    return true;
  }
}

static bool tryDiggingEdge(Adapt* adapter, Collapse& collapse, Entity* e)
{
  Mesh* mesh = adapter->mesh;
  assert(mesh->getType(e) == apf::Mesh::EDGE);
  if ( ! collapse.setEdge(e))
    return false;
  if ( ! collapse.checkClass())
    return false;
  if ( ! collapse.checkTopo())
    return false;
  double q = adapter->input->validQuality;
  if ( ! collapse.tryBothDirections(q))
    return false;
  collapse.destroyOldElements();
  return true;
}

#ifdef DO_FPP
  static bool trySnappingToFPP2(Adapt* a, Collapse& c, Tag* st, Entity* v,
      apf::Up& badElements)
  {
    // make a FPPSnapper Object
    FPPSnapper fpps(a, c, st, v, badElements);
    return fpps.findFPP() ? fpps.snapToFPP() : false;
  }

  static bool trySnappingToFPP(Adapt* a, Collapse& c, Tag* st, Entity* v,
      apf::Up& badElements)
  {
    bool hadItBefore = getFlag(a, v, DONT_COLLAPSE);
    setFlag(a, v, DONT_COLLAPSE);
    bool ok = trySnappingToFPP2(a, c, st, v, badElements);
    if (!hadItBefore)
      clearFlag(a, v, DONT_COLLAPSE);
    return ok;
  }
#endif

static bool tryDigging2(Adapt* a, Collapse& c, apf::Up& badElements)
{
  Mesh* m = a->mesh;
  for (int i = 0; i < badElements.n; ++i) {
    Entity* elem = badElements.e[i];
    Downward edges;
    int nedges;
    nedges = m->getDownward(elem, 1, edges);
    for (int j = 0; j < nedges; ++j){
      /* std::cout << "(i,j) is (" << i << "," << j << ")" << std::endl; */
      if (tryDiggingEdge(a, c, edges[j]))
	return true;
    }
  }
  return false;
}
static bool tryDigging(Adapt* a, Collapse& c, Entity* v,
    apf::Up& badElements)
{
  bool hadItBefore = getFlag(a, v, DONT_COLLAPSE);
  setFlag(a, v, DONT_COLLAPSE);
  bool ok = tryDigging2(a, c, badElements);
  if (!hadItBefore)
    clearFlag(a, v, DONT_COLLAPSE);
  return ok;
}

bool Snapper::run()
{
  dug = false;
  toFPP = false;
  apf::Up badElements;
  bool ok = trySnapping(adapter, snapTag, vert, badElements);
  if (isSimple)
    return ok;
  if (ok)
    return true;
#ifdef DO_FPP
  Mesh* mesh = adapter->mesh;
  assert(mesh->hasTag(vert, snapTag));
  toFPP = trySnappingToFPP(adapter, collapse, snapTag, vert, badElements);
  if (!toFPP)
    return false;
#endif
  dug = tryDigging(adapter, collapse, vert, badElements);
  if (!dug)
    return false;
  return trySnapping(adapter, snapTag, vert, badElements);
}


FPPSnapper::FPPSnapper(Adapt* a, Collapse& c, Tag* st, Entity* v, apf::Up& badElements)
{
  adapter = a;
  snapTag = st;
  collapse = c;
  vert = v;
  problemRegions.n = badElements.n;
  for (int i = 0; i < badElements.n; i++) {
    problemRegions.e[i] = badElements.e[i];
  }
  problemFace = 0;
  problemRegion = 0;
  commEdges.n = 0;
  tol = 1.0e-14;
}

bool FPPSnapper::findFPP()
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

  for (int i = 0; i < n; i++) {
    elem = problemRegions.e[i];
    face = faceOppositeOfVert(elem, vert);
    std::vector<Vector> coords;
    getFaceCoords(face, coords);

    Vector intersect;
    bool isInf;
    bool ok = intersectRayFace(ray, coords, intersect, isInf);
    if (ok){
      if (isInf)
      	std::cout << "Error: could not find a first problem plane!" <<
      	  std::endl;
      Vector newDirection = intersect - ray.start;
      dists.push_back(newDirection.getLength());
      if (dists.back() < minDist) {
      	minDist = dists.back();
      	problemFace = face;
      	problemRegion = elem;
      	intersection = intersect;
      	// do not need to check whether the move is valid since the valid
      	// ones should have been taken care of by this point
      }
      else
      	dists.push_back(minDist + 1.0 + tol);
    }
  }

  //////////
  // this is from old meshAdapt
  // but here since we have digging afterwards there should be
  // no need for this
  //////////
  /* apf::Up coplanarBadElements; */
  /* coplanarBadElements.n = 0; */
  /* if (!problemRegion) { */
  /*   problemRegion = badElements.e[0]; */
  /*   problemFace = faceOppositeOfVert(problemRegion, vert); */
  /*   coplanarBadElements.n = badElements.n; */
  /*   for (int i = 0; i < badElements.n; i++) { */
  /*     coplanarBadElements.e[i] = badElements.e[i]; */
  /*   } */
  /* } */
  /* else */
  /* { */
  /*   minDist += tol; */
  /*   for (int i = 0; i < badElements.n; i++) { */
  /*     if (dists[i] < minDist) { */
  /*     	coplanarBadElements.n++; */
  /*     	coplanarBadElements.e[i] = badElements.e[i]; */
  /*     } */
  /*   } */
  /* } */

  if (!problemRegion)
    return false;

  apf::Up coplanarProblemRegions;
  coplanarProblemRegions.n = 0;
  minDist += tol;
  for (int i = 0; i < problemRegions.n; i++) {
    if (dists[i] < minDist) {
      coplanarProblemRegions.e[coplanarProblemRegions.n] = problemRegions.e[i];
      coplanarProblemRegions.n++;
    }
  }

  findCommonEdges(coplanarProblemRegions);
  return true;
}

bool FPPSnapper::snapToFPP()
{
  bool snapped = false;
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

    double candidateDist = (vCoord - x).getLength();
    if (candidateDist > dist)
      continue;

    // try digging edge here
    snapped = tryDiggingEdge(adapter, collapse, edge);
    if (snapped)
      return true;
  }
  return snapped;
}

Entity* FPPSnapper::faceOppositeOfVert(Entity* e, Entity* v)
{
  Mesh* mesh = adapter->mesh;
  Downward faces;
  Entity* oppositeFace = 0;
  int nDownFaces = mesh->getDownward(e, 2, faces);
  for (int i = 0; i < nDownFaces; i++) {
    Downward verts;
    int nDownVerts = mesh->getDownward(faces[i], 0, verts);
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
  int numDownVerts = mesh->getDownward(oppositeFace, 0, verts);
  bool flag = true;
  for (int i = 0; i < numDownVerts; i++) {
    if (v == verts[i]){
      flag = false;
      break;
    }
  }
  assert(flag);

  return oppositeFace;
}

void FPPSnapper::getFaceCoords(Entity* face, std::vector<Vector>& coords)
{
  Mesh* mesh = adapter->mesh;
  Downward verts;
  int nDownVerts = mesh->getDownward(face, 0, verts);
  assert(nDownVerts);
  for (int i = 0; i < nDownVerts; i++)
    coords.push_back(getPosition(mesh, verts[i]));
}

bool
FPPSnapper::intersectRayFace(const Ray& ray, const std::vector<Vector>& coords,
    Vector& intersection, bool& isInf)
{

  std::cout << ">>>>--------------------" << std::endl;
  /* std::cout << "inside loop at iteration " << i << std::endl; */
  std::cout << "start of ray " << ray.start << std::endl;
  std::cout << "dir   of ray " << ray.dir << std::endl;
  /* std::cout << "number of face verts " << coords.size() << std::endl; */
  /* for (int k = 0; k < (int)coords.size(); k++) { */
  /*   std::cout << "face " << k << ": " << coords[k] << std::endl; */
  /* } */
  std::cout << "--------------------<<<<" << std::endl;


  bool res = false;
  isInf = false;
  if (coords.size() != 3){
    std::cout << "coords.size() is " << coords.size() << std::endl;
    std::cout << "No implementation for non-tri faces!" << std::endl;
    res = false;
  }

  assert(ray.dir.getLength() > tol);
  Vector start = ray.start;
  Vector dir   = ray.dir;

  Vector p0p1 = coords[1] - coords[0];
  Vector p0p2 = coords[2] - coords[0];
  Vector startP0 = coords[0] - start;

  Vector faceAreaVect = apf::cross(p0p1, p0p2);
  double faceAreaSize = faceAreaVect.getLength();

  double vol = dir * faceAreaVect;
  double volPrime = startP0 * faceAreaVect;

  if (vol <= tol * faceAreaSize) {
    if (volPrime <= tol * faceAreaSize) {
      isInf = true;
      res = true;
    }
    else
      res = false;
  }
  else {
    intersection = start - dir * (volPrime / vol);
    Vector newDir = start - intersection;
    if (newDir * dir < 0)
      res = false;
    else
      res = true;
  }
  return res;
}

void FPPSnapper::findCommonEdges(apf::Up& cpRegions)
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
    face = faceOppositeOfVert(region, vert);
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
	if (flag) {
	  commEdges.e[commEdges.n] = edges[i];
	  commEdges.n++;
	}
      }
    }
  }
}

Vector getCenter(Mesh* mesh, Entity* face)
{
  assert(face);
  Downward verts;
  int nDownVerts = mesh->getDownward(face, 0, verts);
  assert(nDownVerts == 3);
  Vector center(0., 0., 0.);
  for (int i = 0; i < nDownVerts; i++)
    center += getPosition(mesh, verts[i]);

  center = center / 3.;
  return center;
}

bool isLowInHigh(Mesh* mesh, Entity* highEnt, Entity* lowEnt)
{
  assert(mesh->getType(highEnt) > mesh->getType(lowEnt));
  Downward down;
  int nDown = mesh->getDownward(highEnt, apf::getDimension(mesh, lowEnt), down);
  for (int i = 0; i < nDown; i++) {
    if (lowEnt == down[i])
      return true;
  }
  return false;
}

}

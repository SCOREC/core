#include "maFirstProblemPlane.h"
#include "lionPrint.h"
#include "pcu_util.h"

namespace ma {

FirstProblemPlane::FirstProblemPlane(Adapt* a, Tag* st)
{
  adapter = a;
  snapTag = st;
  problemFace = 0;
  problemRegion = 0;
  commEdges.clear();
  tol = 1.0e-14;
}

void FirstProblemPlane::setVertex(Entity* v)
{
  vert = v;
}

void FirstProblemPlane::setBadElements(apf::Up& badElements)
{
  for (int i = 0; i < badElements.n; i++) {
    problemRegions.push_back(badElements.e[i]);
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
  size_t n = problemRegions.size();
  Entity* elem;
  Entity* face;
  Ray ray;
  Vector target;
  mesh->getDoubleTag(vert, snapTag, &target[0]);

  ray.start = getPosition(mesh, vert);
  ray.dir   = target - ray.start;

  dists.clear();
  for (size_t i = 0; i < n; i++) {
    elem = problemRegions[i];
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
    problemRegion = problemRegions[0];
    problemFace = getTetFaceOppositeVert(mesh, problemRegion, vert);
    coplanarProblemRegions.n = n;
    for (size_t i = 0; i < n; i++) {
      coplanarProblemRegions.e[i] = problemRegions[i];
    }
  }
  else {
    minDist += tol;
    for (size_t i = 0; i < n; i++) {
      if (dists[i] < minDist) {
        coplanarProblemRegions.e[coplanarProblemRegions.n] = problemRegions[i];
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

  for (size_t i = 0; i < commEdges.size(); i++) {
    edge = commEdges[i];
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
      if (isLowInHigh(mesh, edges[i], vert))
        commEdges.push_back(edges[i]);
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
      if (flag) commEdges.push_back(edges[i]);
    }
  }
}

FirstProblemPlane* getFPP(Adapt* a, Entity* vertex, Tag* snapTag, apf::Up& invalid)
{
  FirstProblemPlane* FPP = new FirstProblemPlane(a, snapTag);
  FPP->setVertex(vertex);
  FPP->setBadElements(invalid);
  std::vector<Entity*> commEdges;
  FPP->getCandidateEdges(commEdges);
  return FPP;
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

//returns the greater index in the case of equality 
static int indexOfMin(double a0, double a1, double a2)
{
  if (a1 < a0) {
    if (a2 < a1) return 2;
    else return 1;
  } else {
    if (a2 < a0) return 2;
    else return 0;
  }
}

static Vector projOnTriPlane(Adapt* a, Entity* vert, Vector normal, Vector v0)
{
  double magN = normal*normal;
  Vector vertPos = getPosition(a->mesh, vert);
  double magCP = (vertPos-v0) * normal;
  double ratio=magCP/magN;

  Vector result;
  for (int i=0; i<3; ++i)
    result[i]=vertPos[i]-ratio*normal[i];
  
  return result;
}

/*
  Given a poorly-shaped tetrahedron, a base triangle and the opposite vertex of the base,
  determine the following information:
  1. the key mesh entities to apply local mesh modification
  2. area of the four face

  return 0 : if an edge is degenerated
     3,5,6 : the tetrahedron has two large dihedral angles. The opposite edges will be stored in ents[0], ents[1].
   1,2,4,7 : the tetrahedron has three large angles. The largeest face is stored in ents[0].
*/
int getTetStats(Adapt* a, Entity* vert, Entity* face, Entity* region, Entity* ents[4], double area[4])
{
  Entity* faceEdges[3];
  a->mesh->getDownward(face, 1, faceEdges);

  Entity* verts[3];
  a->mesh->getDownward(face, 0, verts);

  Vector facePos[3];
  Entity* edges[6];
  for (int i=0; i<3; i++) {
    edges[i]=faceEdges[i];
    facePos[i]=getPosition(a->mesh, verts[i]);
  }

  Entity* faces[4];
  faces[0]=face;

  Entity* problemFaces[4];
  a->mesh->getDownward(region, 2, problemFaces);
  for (int i=0; i<4; i++) {
    if (problemFaces[i] == face ) continue;
    else if (isLowInHigh(a->mesh, problemFaces[i], edges[0])) faces[1] = problemFaces[i];
    else if (isLowInHigh(a->mesh, problemFaces[i], edges[1])) faces[2] = problemFaces[i];
    else if (isLowInHigh(a->mesh, problemFaces[i], edges[2])) faces[3] = problemFaces[i];
  }

  for (int i=1; i<3; i++) {
    Entity* problemEdges[3];
    a->mesh->getDownward(faces[i], 1, problemEdges);
    for (int j=0; j<3; j++) {
      if (problemEdges[j]==edges[i-1] ) continue;
      else if (isLowInHigh(a->mesh, problemEdges[j], verts[0])) edges[3] = problemEdges[j];
      else if (isLowInHigh(a->mesh, problemEdges[j], verts[1])) edges[4] = problemEdges[j];
      else if (isLowInHigh(a->mesh, problemEdges[j], verts[2])) edges[5] = problemEdges[j];
    }
  }

  /* find normal to the plane */
  Vector v01 = facePos[1] - facePos[0];
  Vector v02 = facePos[2] - facePos[0];
  Vector norm = apf::cross(v01, v02);

  Vector projection = projOnTriPlane(a, vert, norm, facePos[0]);
  Vector ri = projection - facePos[0];
  Vector rj = projection - facePos[1];
  Vector rk = projection - facePos[2];

  /* determine which side of the edges does the point R lie.
      First get normal vectors */
  Vector normi = apf::cross(v01, ri);
  Vector normj = apf::cross(facePos[2]-facePos[1], rj);
  Vector normk = apf::cross(facePos[0]-facePos[2], rk);

  Vector mag;
  mag[0]=normi*norm;
  mag[1]=normj*norm;
  mag[2]=normk*norm;

  area[0]=norm*norm;
  area[1]=normi*normi;
  area[2]=normj*normj;
  area[3]=normk*normk;

  int filter[]={1,2,4};
  int bit=0;
  /* examine signs of mag[0], mag[1] and mag[2] */
  for(int i=0; i<3; i++)
    if(mag[i]>0.0)
      bit = bit | filter[i];

  /*  
           010=2   | 011=3  /  001=1
                   |       /
       ------------+--e2--+-----------
                 v0|     /v2
                   | 7  /
           110=6   e0  e1
                   |  /
                   | /    101=5
                   |/
                 v1+
                  /|
                 / |
                  4
  */

 switch( bit ) {
    case 1:{
      int Emap[]={0,4,3};
      int Fmap[]={0,2,3};
      ents[0]=faces[1];
      int i=indexOfMin(area[0],area[2],area[3]);
      ents[1]=edges[Emap[i]];
      ents[2]=faces[Fmap[i]];
      break;
    }
    case 2: {
      int Emap[]={1,4,5};
      int Fmap[]={0,1,3};
      ents[0]=faces[2];
      int i=indexOfMin(area[0],area[1],area[3]);
      ents[1]=edges[Emap[i]];
      ents[2]=faces[Fmap[i]];
      break;   
    }
    case 3: {
      ents[0]=edges[2];
      ents[1]=edges[4];
      ents[2]=((area[0]<area[3]) ? faces[0] : faces[3]);
      ents[3]=((area[1]<area[2]) ? faces[1] : faces[2]);
      break;
    }
    case 4: {
      int Emap[]={2,3,5};
      int Fmap[]={0,1,2};
      ents[0]=faces[3];
      int i=indexOfMin(area[0],area[1],area[2]);
      ents[1]=edges[Emap[i]];
      ents[2]=faces[Fmap[i]];
      break;
    }
    case 5: {
      ents[0]=edges[1];
      ents[1]=edges[3];
      ents[2]=((area[0]<area[2]) ? faces[0] : faces[2]);
      ents[3]=((area[1]<area[3]) ? faces[1] : faces[3]);
      break;
    }
    case 6: {
      ents[0]=edges[0];
      ents[1]=edges[5];
      ents[2]=((area[0]<area[1]) ? faces[0]:faces[1]);
      ents[3]=((area[2]<area[3]) ? faces[2]:faces[3]);
      break;
    }
    case 7: {
      int Emap[]={0,1,2};
      int Fmap[]={1,2,3};
      ents[0]=faces[0];
      int i=indexOfMin(area[1],area[2],area[3]);
      ents[1]=edges[Emap[i]];
      ents[2]=faces[Fmap[i]];
      break;
    }
    default:
      print(a->mesh->getPCU(), "Swap warning: This swap/splt may not work consider more collapses");
  }
  return bit;
}

}
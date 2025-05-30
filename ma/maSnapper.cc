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

Snapper::Snapper(Adapt* a, Tag* st, bool is) : splitCollapse(a), doubleSplitCollapse(a)
{
  adapt = a;
  mesh = a->mesh;
  snapTag = st;
  collapse.Init(a);
  edgeSwap = makeEdgeSwap(a);
  vert = 0;
}

Snapper::~Snapper()
{
  delete edgeSwap;
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
  mesh->getUp(vert,edges);
  apf::Up ovs;
  ovs.n = edges.n;
  for (int i = 0; i < edges.n; ++i)
    ovs.e[i] = apf::getEdgeVertOppositeVert(mesh, edges.e[i], vert);
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

static void flagAndPrint(Adapt* a, Entity* ent, int dim, const char* name)
{
  setFlag(a, ent, CHECKED);
  ma_dbg::dumpMeshWithFlag(a, 0, dim, CHECKED, name, name);
  clearFlag(a, ent, CHECKED);
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

  flagAndPrint(a, FPP->vert, 0, "FPP_Vertex");
  flagAndPrint(a, FPP->problemFace, 2, "FPP_Face");
  flagAndPrint(a, FPP->problemRegion, 3, "FPP_Region");
}

static int indexOfMin(double a0, double a1, double a2)
{
  int k;
  double buf;
  if( a0<a1 )
    { buf=a0; k=0; }
  else
    { buf=a1; k=1; }
  return (buf<a2) ? k:2;
}

static Vector projOnTriPlane(Adapt* a, Entity* face, Entity* vert)
{
  Entity* faceVert[3];
  a->mesh->getDownward(face, 0, faceVert);
  Vector facePos[3];
  for (int i=0; i < 3; ++i)
    facePos[i] = getPosition(a->mesh,faceVert[i]);
  
  Vector normal = apf::cross((facePos[1]-facePos[0]),(facePos[2]-facePos[0]));
  double magN = normal*normal;
  Vector vertPos = getPosition(a->mesh, vert);
  double magCP = (vertPos-facePos[0]) * normal;
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
  2. area of the four faces
  3. the intersection of two intersected opposite edges in case two large dihedral angles
  
  return 0   : if an edge is degenerated
          1-7 : the index indicating the location of projection point
*/
static int getTetStats(Adapt* a, FirstProblemPlane* FPP, Entity* ents[4], double area[4]) 
{
  Entity* faceEdges[3];
  a->mesh->getDownward(FPP->problemFace, 1, faceEdges);

  Entity* verts[3];
  a->mesh->getDownward(FPP->problemFace, 0, verts);

  Vector facePos[3];
  Entity* edges[6];
  for (int i=0; i<3; i++) {
    edges[i]=faceEdges[i];
    facePos[i]=getPosition(a->mesh, verts[i]);
  }

  Entity* faces[4];
  faces[0]=FPP->problemFace;

  Entity* problemFaces[4];
  a->mesh->getDownward(FPP->problemRegion, 2, problemFaces);
  for (int i=0; i<4; i++) {
    if (problemFaces[i] == FPP->problemFace ) continue;
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

  Vector projection = projOnTriPlane(a, FPP->problemFace, FPP->vert);
  //TODO: ERROR if projection = any point on problem face

  /* find normal to the plane */
  Vector v01 = facePos[1] - facePos[0];
  Vector v02 = facePos[2] - facePos[0];
  Vector norm = apf::cross(v01, v02);

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
      // double L1_xyz[2][3], L2_xyz[2][3];
      ents[0]=edges[2];
      ents[1]=edges[4];
      ents[2]=((area[0]<area[3]) ? faces[0] : faces[3]);
      ents[3]=((area[1]<area[2]) ? faces[1] : faces[2]);
      // L1_xyz[0][0]=fxyz[0][0];    L1_xyz[1][0]=fxyz[2][0]; 
      // L1_xyz[0][1]=fxyz[0][1];    L1_xyz[1][1]=fxyz[2][1];
      // L1_xyz[0][2]=fxyz[0][2];    L1_xyz[1][2]=fxyz[2][2];
      // L2_xyz[0][0]=fxyz[1][0];    L2_xyz[1][0]=pxyz[0]; 
      // L2_xyz[0][1]=fxyz[1][1];    L2_xyz[1][1]=pxyz[1];
      // L2_xyz[0][2]=fxyz[1][2];    L2_xyz[1][2]=pxyz[2];
      // fromMeshTools::MT_intLineLine2(L1_xyz,L2_xyz,&i,intXYZ); //TODO: This logic is used for double split collapse
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
      // double L1_xyz[2][3], L2_xyz[2][3];
      ents[0]=edges[1];
      ents[1]=edges[3];
      ents[2]=((area[0]<area[2]) ? faces[0] : faces[2]);
      ents[3]=((area[1]<area[3]) ? faces[1] : faces[3]);
      // L1_xyz[0][0]=fxyz[1][0];    L1_xyz[1][0]=fxyz[2][0]; 
      // L1_xyz[0][1]=fxyz[1][1];    L1_xyz[1][1]=fxyz[2][1];
      // L1_xyz[0][2]=fxyz[1][2];    L1_xyz[1][2]=fxyz[2][2];
      // L2_xyz[0][0]=fxyz[0][0];    L2_xyz[1][0]=pxyz[0]; 
      // L2_xyz[0][1]=fxyz[0][1];    L2_xyz[1][1]=pxyz[1];
      // L2_xyz[0][2]=fxyz[0][2];    L2_xyz[1][2]=pxyz[2];
      // fromMeshTools::MT_intLineLine2(L1_xyz,L2_xyz,&i,intXYZ); //TODO: This logic is used for double split collapse
      break;
    }
    case 6: {
      // double L1_xyz[2][3], L2_xyz[2][3];
      ents[0]=edges[0];
      ents[1]=edges[5];
      ents[2]=((area[0]<area[1]) ? faces[0]:faces[1]);
      ents[3]=((area[2]<area[3]) ? faces[2]:faces[3]);
      // L1_xyz[0][0]=fxyz[0][0];    L1_xyz[1][0]=fxyz[1][0]; 
      // L1_xyz[0][1]=fxyz[0][1];    L1_xyz[1][1]=fxyz[1][1];
      // L1_xyz[0][2]=fxyz[0][2];    L1_xyz[1][2]=fxyz[1][2];
      // L2_xyz[0][0]=fxyz[2][0];    L2_xyz[1][0]=pxyz[0]; 
      // L2_xyz[0][1]=fxyz[2][1];    L2_xyz[1][1]=pxyz[1];
      // L2_xyz[0][2]=fxyz[2][2];    L2_xyz[1][2]=pxyz[2];
      // fromMeshTools::MT_intLineLine2(L1_xyz,L2_xyz,&i,intXYZ); //TODO: This logic is used for double split collapse
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
      PCU_ALWAYS_ASSERT(false); //TODO: ADD ERROR
  }
  return bit;
}

bool Snapper::trySwapOrSplit(FirstProblemPlane* FPP)
{
  if (FPP->commEdges.n < 2) {
    print(mesh->getPCU(), "Swap failed: consider more collapses before swap");
    return false;
  }
  Entity* ents[4];
  double area[4];
  int bit = getTetStats(adapt, FPP, ents, area);

  double min=area[0];
  for(int i=1; i<4; i++) 
    if( area[i]<min ) min=area[i]; 

  if (area[0]==min) { //TODO: implement
    print(mesh->getPCU(), "Swap failed: small base area case not implemented");
    return false;
  }

  // two large dihedral angles -> key problem: two mesh edges
  if (bit==3 || bit==5 || bit==6) {
    // check edge swapping
    for (int i=0; i<2; i++)
      if (edgeSwap->run(ents[i])) { //TODO: Select best
        numSwap++;
        return true;
      }

    // check split+collapse
    for (int i=0; i<2; i++)
      if (splitCollapse.run(ents[i], FPP->vert, 0)) { //TODO: Select best
        numSplitCollapse++;
        return true;
      }
  }
  // three large dihedral angles -> key entity: a mesh face 
  else {

    // check edge swaps
    Entity* edges[3];
    mesh->getDownward(ents[0], 1, edges);
    for (int i=0; i<3; i++) {
      if (edgeSwap->run(edges[i])) { //TODO: Select best
        numSwap++;
        return true;
      }
    }

    //TODO: IMPLEMENT FACE SWAP

    if (splitCollapse.run(ents[1], FPP->vert, 0)) {
      numSplitCollapse++;
      return true;
    }

    // if (doubleSplitCollapse.run(ents, 0)) { //TODO: TEST DOUBLE SPLIT COLLAPSE
    //   numSplitCollapse++;
    //   return true;
    // }
  }

  print(mesh->getPCU(), "Swap failed: face swap and double split collapse not implemented");
  return false;
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
  if (collapse.setEdge(edge) && collapse.checkClass() && collapse.checkTopo()) {
      collapse.computeElementSets();
      if (collapse.tryThisDirectionNoCancel(0) && collapse.edgesGoodSize()) {
        double quality = getWorstQuality(a, collapse.newElements);
        if (quality > best.quality) {
          best.quality = quality;
          best.edge = edge;
          best.keep = keep;
        }
      }
      collapse.cancel();
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

bool Snapper::tryCollapseTetEdges(FirstProblemPlane* FPP)
{
  apf::Up& commEdges = FPP->commEdges;
  BestCollapse best;

  for (int i=0; i<commEdges.n; i++) {
    Entity* vertex[2];
    mesh->getDownward(commEdges.e[i], 0, vertex);
    for (int j=0; j<2; j++)
      getBestQualityCollapse(adapt, commEdges.e[i], vertex[j], collapse, best);
  }

  for (int i=0; i<commEdges.n; i++) {
    Entity* edge = commEdges.e[i];
    Entity* vertexFPP = getEdgeVertOppositeVert(mesh, edge, vert);
    apf::Up adjEdges;
    mesh->getUp(vertexFPP, adjEdges);
    for (int j=0; j<adjEdges.n; j++) {
      Entity* edgeDel = adjEdges.e[j];
      if (edgeDel==edge) continue;
      if (isLowInHigh(mesh, FPP->problemFace, edgeDel)) continue;
      Entity* vertKeep = getEdgeVertOppositeVert(mesh, edgeDel, vertexFPP);
      if (sameSide(adapt, vertKeep, vert, FPP->problemFace)) continue;
      getBestQualityCollapse(adapt, edgeDel, vertKeep, collapse, best);
    }
  }

  if (best.quality > 0) {
    numCollapse++;
    return tryCollapseEdge(adapt, best.edge, best.keep, collapse);
  }
  else return false;
}

bool Snapper::tryReduceCommonEdges(FirstProblemPlane* FPP)
{
  apf::Up& commEdges = FPP->commEdges;
  BestCollapse best;

  Entity* pbEdges[3];
  mesh->getDownward(FPP->problemFace, 1, pbEdges);
  switch(commEdges.n) {
    case 2: {
      Entity* v1 = getEdgeVertOppositeVert(mesh, commEdges.e[0], vert);
      Entity* v2 = getEdgeVertOppositeVert(mesh, commEdges.e[1], vert);
      
      for (int i=0; i<3; i++) {
        Entity* pbVert[2];
        mesh->getDownward(pbEdges[i], 0, pbVert);
        if (pbVert[0] == v1 && pbVert[1] == v2) break;
        if (pbVert[1] == v1 && pbVert[0] == v2) break;
        for (int j=0; j<2; j++)
          getBestQualityCollapse(adapt, pbEdges[i], pbVert[j], collapse, best);
      }
      break;
    }
    case 3: {
      for (int i=0; i<3; i++) {
        Entity* pbVert[2];
        mesh->getDownward(pbEdges[i], 0, pbVert);
        for (int j=0; j<2; j++)
          getBestQualityCollapse(adapt, pbEdges[i], pbVert[j], collapse, best);
      }
      break;
    }
  }
  if (best.quality > 0) {
    numCollapse++;
    return tryCollapseEdge(adapt, best.edge, best.keep, collapse);
  }
  else return false;
}

bool Snapper::tryCollapseToVertex(FirstProblemPlane* FPP)
{
  Vector position = getPosition(mesh, vert);
  Vector target;
  mesh->getDoubleTag(vert, snapTag, &target[0]);
  double distTarget = (position - target).getLength();

  BestCollapse best;

  for (size_t i = 0; i < FPP->commEdges.n; ++i) {
    Entity* edge = FPP->commEdges.e[i];
    Entity* vertexOnFPP = getEdgeVertOppositeVert(mesh, edge, vert);
    Vector vFPPCoord = getPosition(mesh, vertexOnFPP);
    //TODO: add logic for boundary layers
    double distToFPPVert = (vFPPCoord - target).getLength();
    if (distToFPPVert > distTarget) continue;
    getBestQualityCollapse(adapt, edge, vert, collapse, best);
  }

  if (best.quality > 0) {
    numCollapseToVtx++;
    return tryCollapseEdge(adapt, best.edge, best.keep, collapse);
  }
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
  return tryReposition(adapt, vert, snapTag, invalid);
}

bool Snapper::run()
{
  apf::Up invalid;
  bool success = tryReposition(adapt, vert, snapTag, invalid);

  if (success) {
    numSnapped++;
    mesh->removeTag(vert,snapTag);
    clearFlag(adapt, vert, SNAP);
    return true;
  }

  FirstProblemPlane* FPP = getFPP(adapt, vert, snapTag, invalid);

  if (!success) success = tryCollapseToVertex(FPP);
  // if (!success) success = tryReduceCommonEdges(FPP);
  if (!success) success = tryCollapseTetEdges(FPP);
  if (!success) success = trySwapOrSplit(FPP);

  if (!success) {
    numFailed++;
    mesh->removeTag(vert,snapTag);
    clearFlag(adapt, vert, SNAP);
  }
  // if (!success && numFailed == 1) printFPP(adapt, FPP);
  
  if (FPP) delete FPP;
  return success;
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

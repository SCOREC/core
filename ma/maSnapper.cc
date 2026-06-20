/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
/**
 * \file maSnapper.cc
 * \brief Definition of maSnapper.h file.
 * This file contains functions to move a point to the model surface. As described
 * in Li's thesis it will first try to collapse in the target direction. Otherwise 
 * it will collapse to simplify the region and attempt other operators such as
 * swap, split collapse, double split collapse.
*/
/**
 * \file maSnapper.cc
 * \brief Definition of maSnapper.h file.
 * This file contains functions to move a point to the model surface. As described
 * in Li's thesis it will first try to collapse in the target direction. Otherwise 
 * it will collapse to simplify the region and attempt other operators such as
 * swap, split collapse, double split collapse.
*/
#include "maSnapper.h"
#include "maAdapt.h"
#include "maShapeHandler.h"
#include "maFaceSwap.h"
#include "maSnap.h"
#include "maDBG.h"
#include <apfCavityOp.h>
#include <pcu_util.h>
#include <lionPrint.h>
#include <iostream>
#include "apfGeometry.cc"

namespace ma {

Snapper::Snapper(Adapt* a, Tag* st) : mesh(a->mesh), splitCollapse(a), doubleSplitCollapse(a), reposition(a)
{
  adapt = a;
  adapt = a;
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

//Write snapping data to files for debugging purposes
//In order to view relevant information it is neccessary to hide entities with relevent flag in vtk viewer
#if defined(DEBUG_FPP)
static void flagAndPrint(Adapt* a, Entity* ent, int dim, const char* name)
{
  setFlag(a, ent, CHECKED);
  ma_dbg::dumpMeshWithFlag(a, 0, dim, CHECKED, name, name);
  clearFlag(a, ent, CHECKED);
}
#endif

//Write snapping data to files for debugging purposes
#if defined(DEBUG_FPP)
static void printFPP(Adapt* a, FirstProblemPlane* FPP)
{
  ma_dbg::useFieldInfo(a, [a, FPP] {
    apf::writeVtkFiles("FPP_Mesh", a->mesh);
    EntitySet invalid;
    for (auto e : FPP->problemRegions) invalid.insert(e);
    ma_dbg::createCavityMesh(a, invalid, "FPP_Invalid");

    for (auto e : FPP->commEdges) setFlag(a, e, CHECKED);
    ma_dbg::dumpMeshWithFlag(a, 0, 1, CHECKED, "FPP_CommEdges", "FPP_CommEdges");
    for (auto e : FPP->commEdges) clearFlag(a, e, CHECKED);

    flagAndPrint(a, FPP->vert, 0, "FPP_Vertex");
    flagAndPrint(a, FPP->problemFace, 2, "FPP_Face");
    flagAndPrint(a, FPP->problemRegion, 3, "FPP_Region");

    EntitySet adjacent1 = getNextLayer(a, invalid);
    ma_dbg::createCavityMesh(a, adjacent1, "FPP_ADJACENT_1");
    EntitySet adjacent2 = getNextLayer(a, adjacent1);
    ma_dbg::createCavityMesh(a, adjacent2, "FPP_ADJACENT_2");
  });
}
#endif

/*
  We perform this last to make sure that we have a simple region where we can determine the
  best operation to perform and because we want to avoid creating more vertices to snap since
  we could get stuck in an infinite loop of creating and snapping those vertices.
*/
bool Snapper::trySwapOrSplit(FirstProblemPlane* FPP)
{
  Entity* ents[4] = {0};
  double area[4];
  auto problemType = getTetStats(adapt, FPP->vert, FPP->problemFace, FPP->problemRegion, ents, area);
  double qual = adapt->input->validQuality;

  double min=area[0];
  for(int i=1; i<4; i++) 
    if( area[i]<min ) min=area[i]; 

  if (area[0]==min) {
    Entity* edges[3];
    mesh->getDownward(FPP->problemFace, 1, edges);
    Entity* longest = edges[0];
    for (int i=1; i<3; i++)
      if (adapt->sizeField->measure(edges[i]) > adapt->sizeField->measure(longest))
        longest = edges[i];

    if (edgeSwap->run(longest)) { numSwap++; return true; }
    if (splitCollapse.run(longest, FPP->vert, qual)) { numSplitCollapse++; return true; }
  }

  if (ents[0] == 0)
    return false;

  // two large dihedral angles -> key problem: two mesh edges
  if (problemType == ProblemType::TWOLARGEANGLES) {
    if (edgeSwap->run(ents[0])) { numSwap++; return true; }
    if (edgeSwap->run(ents[1])) { numSwap++; return true; }
    if (splitCollapse.run(ents[0], FPP->vert, qual)) { numSplitCollapse++; return true; }
    if (splitCollapse.run(ents[1], FPP->vert, qual)) { numSplitCollapse++; return true; }
    if (doubleSplitCollapse.run(ents, qual)) { numSplitCollapse++; return true; }
  }
  // three large dihedral angles -> key entity: a mesh face
  else {
    Entity* edges[3];
    mesh->getDownward(ents[0], 1, edges);
    for (int i=0; i<3; i++) 
      if (edgeSwap->run(edges[i])) { numSwap++; return true; }
    if (splitCollapse.run(ents[1], FPP->vert, qual)) { numSplitCollapse++; return true; }
  }
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
      collapse.tryBothDirections(a->input->validQuality)) {
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

/*
  Perfomes a collapse operation and stores the operation in best if the quality is better, 
  then cancels the collapse. We want to pick the highest quality after collapsing to the
  first problem plane so future operations are more likely to succeed.
*/
static void getBestQualityCollapse(Adapt* a, Entity* edge, Entity* keep, Collapse& collapse, BestCollapse& best)
{
  PCU_ALWAYS_ASSERT(a->mesh->getType(edge) == apf::Mesh::EDGE);
  bool alreadyFlagged = true;
  if (keep) alreadyFlagged = getFlag(a, keep, DONT_COLLAPSE);
  if (!alreadyFlagged) setFlag(a, keep, DONT_COLLAPSE);
  if (collapse.setEdge(edge) && collapse.checkClass() && collapse.checkTopo()) {
      collapse.computeElementSets();
      if (collapse.tryThisDirectionNoCancel(a->input->validQuality) && collapse.edgesGoodSize()) {
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

//returns if testVert and refVert are on the same side of the face
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

/*
  If collapsing the common edges failed we want to try collapsing any edge that will
  move us towards the first problem plane. We try collapses first in order to simplify
  the region until we can perform smarter operations.
*/
bool Snapper::tryCollapseTetEdges(FirstProblemPlane* FPP)
{
  std::vector<Entity*>& commEdges = FPP->commEdges;
  BestCollapse best;

  for (size_t i=0; i<commEdges.size(); i++) {
    Entity* vertex[2];
    mesh->getDownward(commEdges[i], 0, vertex);
    for (int j=0; j<2; j++)
      getBestQualityCollapse(adapt, commEdges[i], vertex[j], collapse, best);
  }

  for (size_t i=0; i<commEdges.size(); i++) {
    Entity* edge = commEdges[i];
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

/*
  If collapsing to the first problem plane has failed then we want
  to collapse edges on the first problem plane in order to simplify
  the region future operations are more likely to succeed.
*/
bool Snapper::tryReduceCommonEdges(FirstProblemPlane* FPP)
{
  std::vector<Entity*>& commEdges = FPP->commEdges;
  BestCollapse best;

  Entity* pbEdges[3];
  mesh->getDownward(FPP->problemFace, 1, pbEdges);
  switch(commEdges.size()) {
    case 2: {
      Entity* v1 = getEdgeVertOppositeVert(mesh, commEdges[0], vert);
      Entity* v2 = getEdgeVertOppositeVert(mesh, commEdges[1], vert);
      
      for (int i=0; i<3; i++) {
        Entity* pbVert[2];
        mesh->getDownward(pbEdges[i], 0, pbVert);
        if (pbVert[0] == v1 && pbVert[1] == v2) continue;
        if (pbVert[1] == v1 && pbVert[0] == v2) continue;
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

/*
  First we try collapsing to a vertex on the first problem plane because this is the most likely
  operation to succeed. We first try collapsing to the common edges among the invalid regions 
  since those are more likely to succeed.
*/
bool Snapper::tryCollapseToVertex(FirstProblemPlane* FPP)
{
  Vector position = getPosition(mesh, vert);
  Vector target;
  mesh->getDoubleTag(vert, snapTag, &target[0]);
  double distTarget = (position - target).getLength();

  BestCollapse best;

  for (size_t i = 0; i < FPP->commEdges.size(); ++i) {
    Entity* edge = FPP->commEdges[i];
    Entity* vertexOnFPP = getEdgeVertOppositeVert(mesh, edge, vert);
    Vector vFPPCoord = getPosition(mesh, vertexOnFPP);
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

bool Snapper::trySimpleSnap()
{
  if (!mesh->hasTag(vert, snapTag)) return true;
  Vector target;
  adapt->mesh->getDoubleTag(vert, snapTag, &target[0]);
  return reposition.move(vert, target);
}
#if defined(DEBUG_FPP)
static int DEBUGFAILED=0;
#endif

/*
This function will attempt to move vert to the model surface, if it can not do so then it
will atleast move to the first problem plane as described in Li's thesis. It might take multiple
iterations for vert to reach the model surface. Li's thesis was missing some details on how
to apply certain opperators so other algoritms were adapted from old scorec libraries.
*/
bool Snapper::run()
{
  if (!mesh->hasTag(vert, snapTag)) return true;
  Vector target;
  adapt->mesh->getDoubleTag(vert, snapTag, &target[0]);
  bool success = reposition.move(vert, target);
  apf::Up& invalid = reposition.getInvalid();

  if (success) {
    numSnapped++;
    mesh->removeTag(vert,snapTag);
    clearFlag(adapt, vert, SNAP);
    return true;
  }

  FirstProblemPlane* FPP=0;
  if (mesh->getDimension() == 3) {
    if (!success) FPP = getFPP(adapt, vert, target, invalid);
    if (!success) success = tryCollapseToVertex(FPP);
    if (!success) success = tryReduceCommonEdges(FPP);
    if (!success) success = tryCollapseTetEdges(FPP);
    if (!success) success = trySwapOrSplit(FPP);
  }

  if (!success) numFailed++;
  #if defined(DEBUG_FPP)
  if (!success && ++DEBUGFAILED == 1) printFPP(adapt, FPP);
  #endif
  if (FPP) delete FPP;
  return success;
}

EntitySet getNextLayer(Adapt* a, EntitySet& tets)
{
  EntitySet adjacent;
    APF_ITERATE(ma::EntitySet,tets,it) {
      Entity* faces[4];
      a->mesh->getDownward(*it, 2, faces);
      for (int f=0; f<4; f++) {
        apf::Up nextLayer;
        a->mesh->getUp(faces[f], nextLayer);
        for (int n=0; n<nextLayer.n; n++)
          adjacent.insert(nextLayer.e[n]);
      }
    }
    return adjacent;
}

}

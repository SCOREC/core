/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include <iostream>

#include <PCU.h>
#include <lionPrint.h>

#include "em.h"

#include <apfMesh.h>
#include <apfShape.h>
#include <apfField.h>
#include <apfCavityOp.h>
#include "crv.h"
#include "crvShape.h"

#include <mthQR.h>
#include <mth.h>
#include <mth_def.h>

#include <limits>

#include <set>
#include <pcu_util.h>

namespace em {

/* overall information useful during equilibration */
struct Equilibration {
  apf::Mesh* mesh;
  /* mesh dimension, so far handling 3 only */
  int dim;
  /* polynomial order of Nedelec space */
  int order;
  /* input scalar field containing Nedelec dofs for electric field */
  apf::Field* ef;
};

static void setupEquilibration(Equilibration* eq, apf::Field* f)
{
  eq->mesh = apf::getMesh(f);
  eq->dim = eq->mesh->getDimension();
  eq->ef = f;
  eq->order = f->getShape()->getOrder();
}

struct QRDecomp {
  mth::Matrix<double> Q;
  mth::Matrix<double> R;
};

typedef std::set<apf::MeshEntity*> EntitySet;

struct EdgePatch {
  apf::Mesh* mesh;
  Equilibration* equilibration;
  /* the entity around which the patch
     is centered. a patch collects entities
     (faces & tets) around this edge entity */
  apf::MeshEntity* entity;
  bool isOnBdry;
  EntitySet tets;
  EntitySet faces;
  mth::Matrix<double> A;
  mth::Vector<double> b;
  QRDecomp qr;

};

static void setupEdgePatch(EdgePatch* p, Equilibration* eq)
{
  p->mesh = eq->mesh;
  p->equilibration = eq;
  p->entity = 0;
}

static void startEdgePatch(EdgePatch* p, apf::MeshEntity* e)
{
  p->tets.clear();
  p->faces.clear();
  p->entity = e;
  p->isOnBdry = crv::isBoundaryEntity(p->mesh, p->entity);
}

static void addEntityToPatch(EdgePatch* p, apf::MeshEntity* e)
{
  if(p->mesh->getType(e) == apf::Mesh::TRIANGLE)
    p->faces.insert(e);
  if(p->mesh->getType(e) == apf::Mesh::TET)
    p->tets.insert(e);
}

static void addEntitiesToPatch(EdgePatch* p, apf::DynamicArray<apf::MeshEntity*>& es)
{
  for (std::size_t i=0; i < es.getSize(); ++i)
    addEntityToPatch(p, es[i]);
}

static bool getInitialEdgePatch(EdgePatch* p, apf::CavityOp* o)
{
  if ( ! o->requestLocality(&p->entity,1))
    return false;
  apf::DynamicArray<apf::MeshEntity*> adjacent;
  p->mesh->getAdjacent(p->entity, 3, adjacent);
  addEntitiesToPatch(p, adjacent);
  p->mesh->getAdjacent(p->entity, 2, adjacent);
  addEntitiesToPatch(p, adjacent);
  return true;
}

static bool buildEdgePatch(EdgePatch* p, apf::CavityOp* o)
{
  if (!getInitialEdgePatch(p, o)) return false;
  return true;
}

static void assembleLHS(EdgePatch* p)
{
  int ne = p->tets.size();
  int nf = p->faces.size();
  printf("ne %d nf %d \n", ne, nf); // TODO remove
  if( crv::isBoundaryEntity(p->mesh, p->entity) ) {
    p->A.resize(ne+nf, ne+nf);
    p->A.zero();
    for (int i = 0; i < nf; i++)
      p->A(i,i) = 2.;
    for (int i = 0; i < ne-1; i++) {
      p->A(i+nf,i) = 1.; p->A(i+nf,i+1) = -1.;
      p->A(i,i+nf) = 1.; p->A(i+1,i+nf) = -1.;
    }
    p->A(ne+nf-1, ne-1) = 1.; p->A(ne+nf-1, ne) = 1.;
    p->A(ne-1, ne+nf-1) = 1.; p->A(ne, ne+nf-1) = 1.;

    std::cout << "boundary" << std::endl; // TODO remove
    std::cout << p->A << std::endl; // TODO remove
  }
  else if( ! crv::isBoundaryEntity(p->mesh, p->entity) ) {
    mth::Matrix<double> m(ne, nf);
    m.zero();
    for (int i = 0; i < ne-1; i++) {
      m(i,i) = 1.;  m(i,i+1) = -1.;
    }
    m(ne-1,0) = -1.; m(ne-1,ne-1) = 1.;

    // m is singular so do (m*mt) + 1.0
    // to pick a particular solution
    mth::Matrix<double> mt(nf,ne);
    mth::transpose(m, mt);
    mth::multiply(m, mt, p->A);

    for (int i = 0; i < ne; i++)
      for (int j = 0; j < ne; j++)
        p->A(i,j) += 1.;
    std::cout << "interior" << std::endl; // TODO remove
    std::cout << p->A << std::endl; // TODO remove
  }
}

// computes local bilinear form integral restricted to
// an edge of an element
static int getLocalBLFIntegral(EdgePatch* p, apf::MeshEntity* tet)
{
  // find position of edge in downward edges of element
  apf::Downward e;
  int ne = p->mesh->getDownward(tet, 1, e);
  int ei = apf::findIn(e, ne, p->entity);
  // get Element Dofs
  //
  // assemble curl curl element matrix
  // assemble vector mass element matrix
  // add element matrices
  // multiply element matrix with element dofs
  // pick edge index from the resulting vector
}


static void assembleRHS(EdgePatch* p)
{
  if (p->isOnBdry) {
    p->b.resize(p->tets.size() + p->faces.size());
    p->b.zero();
  }
  else {
    p->b.resize(p->tets.size());
    p->b.zero();
  }

  apf::MeshEntity* edge = p->entity[i];
  int ne = p->tets.size();
  for (int i = 0; i < ne; i++) {
    apf::MeshEntity* tet = p->tets[i];
    int blfIntegral = getLocalBLFIntegral(p, tet);
  }
  // TODO computeBilinearFormIntegral
  // TODO computeLinearFormIntegral
  // TODO computeFluxTermIntegral


}

// The following two functions help order tets and faces in a cavity in a
// clockwise (or ccw) direction.
static apf::MeshEntity* getTetOppFaceSharingEdge(
    apf::Mesh* m, apf::MeshEntity* t, apf::MeshEntity* f, apf::MeshEntity* e)
{
  apf::MeshEntity* fs[4];
  m->getDownward(t, 2, fs);
  for (int i = 0; i < 4; i++) {
    if (fs[i] == f) continue;
    apf::MeshEntity* es[3];
    m->getDownward(fs[i], 1, es);
    if (apf::findIn(es, 3, e) > -1)
      return fs[i];
  }
  return 0;
}
static void getOrderedTetsandFaces(apf::Mesh* mesh, apf::MeshEntity* edge,
     EntitySet& tets, EntitySet& faces)
{
  tets.clear();
  faces.clear();
  if( ! crv::isBoundaryEntity(mesh, edge) ) { // interior edge patches
    apf::MeshEntity* currentFace = mesh->getUpward(edge, 0);
    apf::Up up;
    mesh->getUp(currentFace, up);
    PCU_ALWAYS_ASSERT(up.n == 2);
    apf::MeshEntity* firstTet = up.e[0];
    apf::MeshEntity* nextTet  = up.e[1];
    tets.insert(firstTet);
    apf::MeshEntity* firstFace = getTetOppFaceSharingEdge(mesh, firstTet,
        currentFace, edge);
    faces.insert(firstFace);

    while (nextTet != firstTet) {
      tets.insert(nextTet);
      faces.insert(currentFace);
      currentFace = getTetOppFaceSharingEdge(mesh, nextTet, currentFace, edge);
      PCU_ALWAYS_ASSERT(currentFace);
      apf::Up up;
      mesh->getUp(currentFace, up);
      PCU_ALWAYS_ASSERT(up.n == 2);
      if (nextTet != up.e[0])
        nextTet = up.e[0];
      else
        nextTet = up.e[1];
    }
  }
  else { // boundary edge patches
    apf::Up up;
    mesh->getUp(edge, up);
    apf::MeshEntity* firstFace;
    for (int i = 0; i < up.n; i++) {
      if ( crv::isBoundaryEntity(mesh, up.e[i]) ) {
        firstFace = up.e[i]; break;
      }
    }
    faces.insert(firstFace);
    mesh->getUp(firstFace, up);
    PCU_ALWAYS_ASSERT(up.n == 1);
    apf::MeshEntity* firstTet = up.e[0];
    tets.insert(firstTet);

    apf::MeshEntity* nextFace = getTetOppFaceSharingEdge(mesh, firstTet,
        firstFace, edge);
    apf::MeshEntity* nextTet = firstTet;
    mesh->getUp(nextFace, up);
    while( up.n == 2) {
      faces.insert(nextFace);
      if (nextTet != up.e[0])
        nextTet = up.e[0];
      else
        nextTet = up.e[1];
      tets.insert(nextTet);

      nextFace = getTetOppFaceSharingEdge(mesh, nextTet, nextFace, edge);
      mesh->getUp(nextFace, up);
    }
    faces.insert(nextFace);
  }
}

static void runErm(EdgePatch* p)
{
  assembleLHS(p);
  // TODO decompose A into Q and R
  //std::vector<apf::MeshEntity*> otets; // ordered tets // TODO remove these two
  //std::vector<apf::MeshEntity*> ofaces; // ordered faces
  getOrderedTetsandFaces(p->mesh, p->entity, p->tets, p->faces); // TODO maybe use p->tets, p->faces


  // TODO assemble RHS
}


class EdgePatchOp : public apf::CavityOp
{
public:
  EdgePatchOp(Equilibration* eq):
    apf::CavityOp(eq->mesh)
  {
    setupEdgePatch(&edgePatch, eq);
  }
  virtual Outcome setEntity(apf::MeshEntity* e)
  {
    startEdgePatch(&edgePatch, e);
    if ( ! buildEdgePatch(&edgePatch, this))
      return REQUEST;
    return OK;
    return SKIP;
  }
  virtual void apply()
  {
    runErm(&edgePatch);
  }
  EdgePatch edgePatch;
};

void equilibrateResiduals(apf::Field* f)
{
  Equilibration equilibration;
  setupEquilibration(&equilibration, f);
  EdgePatchOp op(&equilibration);
  op.applyToDimension(1); // edges
}







}

/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include <iostream>

#include <apf.h>
#include <PCU.h>
#include <lionPrint.h>

#include "em.h"

#include <apfMesh.h>
#include <apfShape.h>
#include <apfField.h>
#include <apfElement.h>
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

static void assembleCurlCurlElementMatrix(apf::Mesh* mesh,apf::MeshEntity* e,
    apf::Field* f, mth::Matrix<double>& elmat)
{
  apf::FieldShape* fs = f->getShape();
  int type = mesh->getType(e);
  int nd = apf::countElementNodes(fs, type);
  int dim = apf::getDimension(mesh, e);
  int dimc = (dim == 3) ? 3 : 1;
  double w;

  apf::NewArray<apf::Vector3> curlshape(nd);
  mth::Matrix<double> phys_curlshape(nd, dimc);
  elmat.resize(nd,nd);

  apf::MeshElement* me = apf::createMeshElement(mesh, e);
  apf::Element* el = apf::createElement(f, me);
  int order = 2 * fs->getOrder() - 2;
  int np = apf::countIntPoints(me, order); // int points required

  elmat.zero();
  apf::Vector3 p;
  for (int i = 0; i < np; i++) {
    apf::getIntPoint(me, order, i, p);
    double weight = apf::getIntWeight(me, order, i);
    apf::Matrix3x3 J;
    apf::getJacobian(me, p, J);
    double jdet = apf::getJacobianDeterminant(J, dim);
    w = weight / jdet;

    if (dim == 3) {
      el->getShape()->getLocalVectorCurls(mesh, e, p, curlshape);
      phys_curlshape.zero();
      for (int i = 0; i < nd; i++)
        for (int j = 0; j < dim; j++)
          for (int k = 0; k < dim; k++)
            phys_curlshape(i,j) += curlshape[i][k] * J[k][j];
    }
    else {
      /*
      TODO 2D TODO in apfNedelec.cc
      el->getShape()->getLocalVectorCurls(mesh, e, p, curlshape);
      phys_curlshape.zero();
      for (int i = 0; i < nd; i++)
        for (int j = 0; j < dimc; j++)
        */
    }
    mth::Matrix<double> phys_curlshapeT;
    mth::transpose(phys_curlshape, phys_curlshapeT);
    mth::Matrix<double> M (nd, nd);
    M.zero();
    mth::multiply(phys_curlshape, phys_curlshapeT, M);
    M *= w;
    elmat += M;
  }
}

static void assembleVectorMassElementMatrix(apf::Mesh* mesh,apf::MeshEntity* e,
    apf::Field* f, mth::Matrix<double>& elmat)
{
  apf::FieldShape* fs = f->getShape();
  int type = mesh->getType(e);
  int nd = apf::countElementNodes(fs, type);
  int dim = apf::getDimension(mesh, e);
  double w;

  apf::NewArray<apf::Vector3> vectorshape(nd);
  elmat.resize(nd,nd);

  apf::MeshElement* me = apf::createMeshElement(mesh, e);
  apf::Element* el = apf::createElement(f, me);
  int order = 2 * fs->getOrder();
  int np = apf::countIntPoints(me, order); // int points required

  elmat.zero();
  apf::Vector3 p;
  for (int i = 0; i < np; i++) {
    apf::getIntPoint(me, order, i, p);
    double weight = apf::getIntWeight(me, order, i);
    apf::Matrix3x3 J;
    apf::getJacobian(me, p, J);
    double jdet = apf::getJacobianDeterminant(J, dim);
    w = weight * jdet;

    apf::getVectorShapeValues(el, p, vectorshape);
    mth::Matrix<double> vectorShape (nd, dim);
    for (int i = 0; i < nd; i++)
      for (int j = 0; j < dim; j++)
        vectorShape(i,j) = vectorshape[i][j];

    mth::Matrix<double> vectorShapeT (dim, nd);
    mth::transpose(vectorShape, vectorShapeT);
    mth::Matrix<double> M (nd,nd);
    M.zero();
    mth::multiply(vectorShape, vectorShapeT, M);
    M *= w;
    elmat += M;
  }
}

// computes local bilinear form integral restricted to
// an edge of an element
static double getLocalBLFIntegral(EdgePatch* p, apf::MeshEntity* tet)
{
  // find position of edge in downward edges of element
  apf::Downward e;
  int ne = p->mesh->getDownward(tet, 1, e);
  int ei = apf::findIn(e, ne, p->entity);
  // get Element Dofs
  apf::MeshElement* me = apf::createMeshElement(p->mesh, tet);
  apf::Element* el = apf::createElement(p->equilibration->ef, me);
  int type = p->mesh->getType(tet);
  int nd = apf::countElementNodes(el->getFieldShape(), type);
  apf::NewArray<double> d (nd);
  el->getElementDofs(d);
  mth::Vector<double> dofs (nd);
  for (int i = 0; i < nd; i++) // TODO cleanup
    dofs(i) = d[i];
  // assemble curl curl element matrix
  mth::Matrix<double> curl_elmat;
  assembleCurlCurlElementMatrix(p->mesh, tet,
      p->equilibration->ef, curl_elmat);
  // assemble vector mass element matrix
  mth::Matrix<double> mass_elmat;
  assembleVectorMassElementMatrix(p->mesh, tet,
      p->equilibration->ef, mass_elmat);
  // add element matrices
  mth::Matrix<double> elmat(nd, nd);
  elmat.zero();
  elmat += curl_elmat;
  elmat += mass_elmat;
  // multiply element matrix with element dofs
  mth::Vector<double> integrals (nd);
  mth::multiply(elmat, dofs, integrals);
  // pick edge index from the resulting vector
  return integrals(ei);
}

// TODO redo this to allow user access from outside
static void pumiUserFunction(const apf::Vector3& x, mth::Vector<double>& f,
    apf::MeshEntity* tet, apf::Mesh* mesh)
{
  double freq = 1.;
  double kappa = freq * M_PI;
  int dim = apf::getDimension(mesh, tet);
  if (dim == 3) {
      f(0) = (1. + kappa * kappa) * sin(kappa * x[1]);
      f(1) = (1. + kappa * kappa) * sin(kappa * x[2]);
      f(2) = (1. + kappa * kappa) * sin(kappa * x[0]);
  }
  else {
     f(0) = (1. + kappa * kappa) * sin(kappa * x[1]);
     f(1) = (1. + kappa * kappa) * sin(kappa * x[0]);
     f(2) = 0.0;
  }
}

void assembleDomainLFElementVector(apf::Mesh* mesh, apf::MeshEntity* e,
    apf::Field* f, mth::Vector<double>& elvect)
{
  apf::MeshElement* me = apf::createMeshElement(mesh, e);
  apf::Element* el = apf::createElement(f, me);
  int type = mesh->getType(e);
  int nd = apf::countElementNodes(el->getFieldShape(), type);
  int dim = apf::getDimension(mesh, e);

  apf::NewArray<apf::Vector3> vectorshape(nd);
  elvect.resize(nd);
  mth::Vector<double> val (3);
  val.zero();
  apf::Vector3 p;
  double w;

  apf::FieldShape* fs = f->getShape();
  int order = 2 * fs->getOrder();
  int np = apf::countIntPoints(me, order); // int points required

  elvect.zero();
  for (int i = 0; i < np; i++) {
    apf::getIntPoint(me, order, i, p);
    double weight = apf::getIntWeight(me, order, i);
    apf::Matrix3x3 J;
    apf::getJacobian(me, p, J);
    double jdet = apf::getJacobianDeterminant(J, dim);
    w = weight * jdet;

    apf::getVectorShapeValues(el, p, vectorshape);
    pumiUserFunction(p, val, e, apf::getMesh(f)); // eval vector function
    val *= w;

    mth::Matrix<double> vectorShape (nd, dim);
    for (int i = 0; i < nd; i++)
      for (int j = 0; j < dim; j++)
        vectorShape(i,j) = vectorshape[i][j];
    mth::Vector<double> temp (nd);
    temp.zero();
    mth::multiply(vectorShape, val, temp);
    elvect += temp;
  }
}

static double getLocalLFIntegral(EdgePatch* p, apf::MeshEntity* tet)
{
  // find position of edge in downward edges of element
  apf::Downward e;
  int ne = p->mesh->getDownward(tet, 1, e);
  int ei = apf::findIn(e, ne, p->entity);
  // assemble Domain LF Vector
  mth::Vector<double> elvect;
  assembleDomainLFElementVector(p->mesh, tet,
      p->equilibration->ef, elvect);
  return elvect(ei);
}


static void assembleRHS(EdgePatch* p)
{
  if (p->isOnBdry) {
    std::cout << "boundary" << std::endl; // TODO remove
    p->b.resize(p->tets.size() + p->faces.size());
    p->b.zero();
  }
  else {
    std::cout << "interior" << std::endl; // TODO remove
    p->b.resize(p->tets.size());
    p->b.zero();
  }
  printf("ASSEMBLE RHS\n");
  //apf::MeshEntity* edge = p->entity;
  int ne = p->tets.size();
  int nf = p->faces.size();
  for (int i = 0; i < ne; i++) {
    // TODO computeBilinearFormIntegral
    EntitySet::iterator it = std::next(p->tets.begin(), i);
    apf::MeshEntity* tet = *it;
    double blfIntegral = getLocalBLFIntegral(p, tet);
    double lfIntegral = getLocalLFIntegral(p, tet);
    // TODO computeLinearFormIntegral
    // TODO computeFluxTermIntegral
    if(p->isOnBdry)
      p->b(nf+i) = blfIntegral - lfIntegral; // TODO add -flux term integral
    else
      p->b(i) = blfIntegral - lfIntegral;    // TODO add -flux term integral
  }
  std::cout << p->b << std::endl; // TODO remove
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
  assembleRHS(p);

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

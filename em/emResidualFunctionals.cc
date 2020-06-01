/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include <iostream>
#include <cstdlib>
#include <algorithm>

#include <apfCavityOp.h>
#include "apfElement.h"
#include "crv.h"
#include "crvShape.h"

#include "em.h"

using namespace std;

namespace em {

/* overall information useful during equilibration */
struct Equilibration {
  apf::Mesh* mesh;
  /* mesh dimension, so far handling 3 only */
  int dim;
  /* polynomial order of Nedelec space */
  int order;
  /* input scalar field containing Nedelec dofs for solution electric field */
  apf::Field* ef;
  /* output field containing correction values.
   * currently 3 scalar values are stored on each face
   * in order corresponding to downward edges of the face */
  apf::Field* g;
};

static void setupEquilibration(Equilibration* eq, apf::Field* f, apf::Field* g)
{
  eq->mesh = apf::getMesh(f);
  eq->dim = eq->mesh->getDimension();
  eq->ef = f;
  eq->order = f->getShape()->getOrder();
  eq->g = g;

  double zeros[3] = {0., 0., 0.};
  apf::MeshEntity* face;
  apf::MeshIterator* it = eq->mesh->begin(2);
  while ((face = eq->mesh->iterate(it))) {
    apf::setComponents(eq->g, face, 0, zeros);
  }
  eq->mesh->end(it);
}

struct QRDecomp {
  mth::Matrix<double> Q;
  mth::Matrix<double> R;
};

typedef std::vector<apf::MeshEntity*> EntityVector;

struct EdgePatch {
  apf::Mesh* mesh;
  Equilibration* equilibration;
  /* the edge entity around which the patch
     is centered. a patch collects entities
     (faces & tets) around this edge entity */
  apf::MeshEntity* entity;
  bool isOnBdry;
  EntityVector tets;
  EntityVector faces;
  mth::Matrix<double> A;
  mth::Matrix<double> At;
  mth::Matrix<double> T; // T = A*At + 1
  mth::Vector<double> b;
  mth::Vector<double> x;
  QRDecomp qr;
};

static void setupEdgePatch(EdgePatch* ep, Equilibration* eq)
{
  ep->mesh = eq->mesh;
  ep->equilibration = eq;
  ep->entity = 0;
}

static void startEdgePatch(EdgePatch* ep, apf::MeshEntity* e)
{
  ep->tets.clear();
  ep->faces.clear();
  ep->entity = e;
  ep->isOnBdry = crv::isBoundaryEntity(ep->mesh, ep->entity);
}

static void addEntityToPatch(EdgePatch* ep, apf::MeshEntity* e)
{
  if(ep->mesh->getType(e) == apf::Mesh::TRIANGLE)
    ep->faces.push_back(e);
  if(ep->mesh->getType(e) == apf::Mesh::TET)
    ep->tets.push_back(e);
}

static void addEntitiesToPatch(
    EdgePatch* ep, apf::DynamicArray<apf::MeshEntity*>& es)
{
  for (std::size_t i=0; i < es.getSize(); ++i)
    addEntityToPatch(ep, es[i]);
}

static bool getInitialEdgePatch(EdgePatch* ep, apf::CavityOp* o)
{
  if ( ! o->requestLocality(&ep->entity,1))
    return false;
  apf::DynamicArray<apf::MeshEntity*> adjacent;
  ep->mesh->getAdjacent(ep->entity, 3, adjacent);
  addEntitiesToPatch(ep, adjacent);

  ep->mesh->getAdjacent(ep->entity, 2, adjacent);
  addEntitiesToPatch(ep, adjacent);

  return true;
}

static bool buildEdgePatch(EdgePatch* ep, apf::CavityOp* o)
{
  if (!getInitialEdgePatch(ep, o)) return false;
  return true;
}

/*
 * Reference: Leszek Demkowicz, Computing with hp-adaptive
 * finite elements, vol2, equation (11.118, 11.119, 11.122).
 */
static void assembleEdgePatchLHS(EdgePatch* ep)
{
  int ne = ep->tets.size();
  int nf = ep->faces.size();
  if( crv::isBoundaryEntity(ep->mesh, ep->entity) ) {
    ep->T.resize(ne+nf, ne+nf);
    ep->T.zero();
    for (int i = 0; i < nf; i++)
      ep->T(i,i) = 2.;
    for (int i = 0; i < ne-1; i++) {
      ep->T(i+nf,i) = 1.; ep->T(i+nf,i+1) = -1.;
      ep->T(i,i+nf) = 1.; ep->T(i+1,i+nf) = -1.;
    }
    ep->T(ne+nf-1, ne-1) = 1.; ep->T(ne+nf-1, ne) = 1.;
    ep->T(ne-1, ne+nf-1) = 1.; ep->T(ne, ne+nf-1) = 1.;
  }
  else if( ! crv::isBoundaryEntity(ep->mesh, ep->entity) ) {
    ep->A.resize(ne, nf);
    ep->A.zero();
    for (int i = 0; i < ne-1; i++) {
      ep->A(i,i) = 1.;  ep->A(i,i+1) = -1.;
    }
    ep->A(ne-1,0) = -1.; ep->A(ne-1,ne-1) = 1.;
    // A is singular so do (A*At) + 1.0
    // to pick a particular solution
    ep->At.resize(nf,ne);
    mth::transpose(ep->A, ep->At);
    mth::multiply(ep->A, ep->At, ep->T);

    for (int i = 0; i < ne; i++)
      for (int j = 0; j < ne; j++)
        ep->T(i,j) += 1.;
  }
  mth::decomposeQR(ep->T, ep->qr.Q, ep->qr.R);
}

/*
 * Performs Curl Curl integration using curl vector Nedelec shapes
 * TODO Add 2D curl curl integration
 * TODO CLEANUP.
 */
void assembleCurlCurlElementMatrix(apf::Mesh* mesh, apf::MeshEntity* e,
    apf::Field* f, mth::Matrix<double>& elmat)
{
  apf::FieldShape* fs = f->getShape();
  int type = mesh->getType(e);
  PCU_ALWAYS_ASSERT(type == apf::Mesh::TET); // TODO add Triangle
  int nd = apf::countElementNodes(fs, type);
  int dim = apf::getDimension(mesh, e);
  int dimc = (dim == 3) ? 3 : 1;
  double w;

  apf::NewArray<apf::Vector3> curlshape(nd);
  mth::Matrix<double> phys_curlshape(nd, dimc); // TODO remove once curl Piola is in place in apfNedelec
  elmat.resize(nd,nd);

  apf::MeshElement* me = apf::createMeshElement(mesh, e);
  apf::Element* el = apf::createElement(f, me);
  int int_order = 2 * fs->getOrder() - 2;
  int np = apf::countIntPoints(me, int_order);

  elmat.zero();
  apf::Vector3 p;
  for (int i = 0; i < np; i++) {
    apf::getIntPoint(me, int_order, i, p);
    double weight = apf::getIntWeight(me, int_order, i);
    apf::Matrix3x3 J;
    apf::getJacobian(me, p, J);
    double jdet = apf::getJacobianDeterminant(J, dim);
    w = weight / jdet;

    if (dim == 3) { // TODO this is Piola transformation. Put it in apfNedelec
      el->getShape()->getLocalVectorCurls(mesh, e, p, curlshape);
      phys_curlshape.zero();
      for (int j = 0; j < nd; j++)
        for (int k = 0; k < dim; k++)
          for (int l = 0; l < dim; l++)
            phys_curlshape(j,k) += curlshape[j][l] * J[l][k];
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
  apf::destroyElement(el);
  apf::destroyMeshElement(me);
}

/*
 * Performs Vector Vector Mass integration using vector Nedelec shapes
 */
void assembleVectorMassElementMatrix(apf::Mesh* mesh, apf::MeshEntity* e,
    apf::Field* f, mth::Matrix<double>& elmat)
{
  apf::FieldShape* fs = f->getShape();
  int type = mesh->getType(e);
  PCU_ALWAYS_ASSERT(type == apf::Mesh::TET || type == apf::Mesh::TRIANGLE);
  int nd = apf::countElementNodes(fs, type);
  int dim = apf::getDimension(mesh, e);
  int sdim = mesh->getDimension();
  double w;

  apf::NewArray<apf::Vector3> vectorshapes(nd);
  elmat.resize(nd,nd);

  apf::MeshElement* me = apf::createMeshElement(mesh, e);
  apf::Element* el = apf::createElement(f, me);
  int int_order = 2 * fs->getOrder();
  int np = apf::countIntPoints(me, int_order);

  elmat.zero();
  apf::Vector3 p;
  for (int i = 0; i < np; i++) {
    apf::getIntPoint(me, int_order, i, p);
    double weight = apf::getIntWeight(me, int_order, i);
    apf::Matrix3x3 J;
    apf::getJacobian(me, p, J);
    double jdet = apf::getJacobianDeterminant(J, dim);
    w = weight * jdet;

    apf::getVectorShapeValues(el, p, vectorshapes);
    mth::Matrix<double> vectorShapes (nd, sdim);
    for (int j = 0; j < nd; j++)
      for (int k = 0; k < sdim; k++)
        vectorShapes(j,k) = vectorshapes[j][k];

    mth::Matrix<double> vectorShapesT (sdim, nd);
    mth::transpose(vectorShapes, vectorShapesT);

    mth::Matrix<double> M (nd,nd);
    M.zero();
    mth::multiply(vectorShapes, vectorShapesT, M);
    M *= w;
    elmat += M;
  }

  apf::destroyElement(el);
  apf::destroyMeshElement(me);
}

void assembleElementMatrix(apf::Mesh* mesh, apf::MeshEntity*e,
    apf::Field* f, mth::Matrix<double>& elmat)
{
  mth::Matrix<double> curl_elmat, mass_elmat;
  assembleCurlCurlElementMatrix(mesh, e, f, curl_elmat);
  assembleVectorMassElementMatrix(mesh, e, f, mass_elmat);

  elmat.resize(curl_elmat.rows(), curl_elmat.cols());
  elmat.zero();
  elmat += curl_elmat;
  elmat += mass_elmat;
}

/*
 * computes local bilinear form integral restricted to
 * an edge of a tet element
 */
static double getLocalEdgeBLF(EdgePatch* ep, apf::MeshEntity* tet)
{
  PCU_ALWAYS_ASSERT(ep->mesh->getType(tet) == apf::Mesh::TET);
  // findIn edge in downward edges of tet
  apf::Downward e;
  int ne = ep->mesh->getDownward(tet, 1, e);
  int ei = apf::findIn(e, ne, ep->entity);
  // get Element Dofs
  apf::MeshElement* me = apf::createMeshElement(ep->mesh, tet);
  apf::Element* el = apf::createElement(ep->equilibration->ef, me);
  int type = ep->mesh->getType(tet);
  int nd = apf::countElementNodes(el->getFieldShape(), type);
  apf::NewArray<double> d (nd);
  el->getElementDofs(d);
  mth::Vector<double> dofs (nd);
  for (int i = 0; i < nd; i++)
    dofs(i) = d[i];
  // assemble curl curl element matrix
  mth::Matrix<double> curl_elmat;
  assembleCurlCurlElementMatrix(ep->mesh, tet,
      ep->equilibration->ef, curl_elmat);
  // assemble vector mass element matrix
  mth::Matrix<double> mass_elmat;
  assembleVectorMassElementMatrix(ep->mesh, tet,
      ep->equilibration->ef, mass_elmat);
  // add element matrices
  mth::Matrix<double> elmat(nd, nd);
  elmat.zero();
  elmat += curl_elmat;
  elmat += mass_elmat;
  // multiply element matrix with element dofs
  mth::Vector<double> blf_integrals (nd);
  mth::multiply(elmat, dofs, blf_integrals);

  apf::destroyElement(el);
  apf::destroyMeshElement(me);

  // pick edge index from the resulting vector
  // negation of negative ND dofs
  int which, rotate; bool flip;
  apf::getAlignment(ep->mesh, tet, ep->entity, which, flip, rotate);
  if (flip)
    blf_integrals(ei) = -1*blf_integrals(ei);
  return blf_integrals(ei);
}


// TODO QUESTION redo this to allow user access from outside
void pumiUserFunction(apf::Mesh* mesh, apf::MeshEntity* e,
    const apf::Vector3& x, mth::Vector<double>& f)
{
  double freq = 1.;
  double kappa = freq * M_PI;
  int dim = apf::getDimension(mesh, e);
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
  apf::FieldShape* fs = f->getShape();
  int type = mesh->getType(e);
  PCU_ALWAYS_ASSERT(type == apf::Mesh::TET);
  int nd = apf::countElementNodes(fs, type);
  int dim = apf::getDimension(mesh, e);
  double w;

  apf::NewArray<apf::Vector3> vectorshapes(nd);
  elvect.resize(nd);
  mth::Vector<double> val (dim);
  val.zero();
  apf::Vector3 p;

  apf::MeshElement* me = apf::createMeshElement(mesh, e);
  apf::Element* el = apf::createElement(f, me);
  int int_order = 2 * fs->getOrder();
  int np = apf::countIntPoints(me, int_order);

  elvect.zero();
  for (int i = 0; i < np; i++) {
    apf::getIntPoint(me, int_order, i, p);
    double weight = apf::getIntWeight(me, int_order, i);
    apf::Matrix3x3 J;
    apf::getJacobian(me, p, J);
    double jdet = apf::getJacobianDeterminant(J, dim);
    w = weight * jdet;

    apf::getVectorShapeValues(el, p, vectorshapes);
    apf::Vector3 global;
    apf::mapLocalToGlobal(me, p, global);
    pumiUserFunction(mesh, e, global, val);
    val *= w;

    mth::Matrix<double> vectorShapes (nd, dim);
    for (int j = 0; j < nd; j++)
      for (int k = 0; k < dim; k++)
        vectorShapes(j,k) = vectorshapes[j][k];
    mth::Vector<double> V (nd);
    V.zero();
    mth::multiply(vectorShapes, val, V);
    elvect += V;
  }

  apf::destroyElement(el);
  apf::destroyMeshElement(me);
}

/*
 * computes local linear form integral restricted to
 * an edge of a tet element
 */
static double getLocalEdgeLF(EdgePatch* ep, apf::MeshEntity* tet)
{
  PCU_ALWAYS_ASSERT(ep->mesh->getType(tet) == apf::Mesh::TET);
  // findIn edge in downward edges of tet
  apf::Downward e;
  int ne = ep->mesh->getDownward(tet, 1, e);
  int ei = apf::findIn(e, ne, ep->entity);
  // assemble Domain LF Vector
  mth::Vector<double> elvect;
  assembleDomainLFElementVector(ep->mesh, tet,
      ep->equilibration->ef, elvect);
  int which, rotate; bool flip;
  apf::getAlignment(ep->mesh, tet, ep->entity, which, flip, rotate);
  if (flip)
    elvect(ei) = -1*elvect(ei);
  return elvect(ei);
}

/*
 * Given a tet and one of its faces, the vertex of the tet
 * opposite to the given face is returned.
 */
static apf::MeshEntity* getTetOppVert(
    apf::Mesh* m, apf::MeshEntity* t, apf::MeshEntity* f)
{
  apf::Downward fvs;
  int fnv = m->getDownward(f, 0, fvs);
  apf::Downward tvs;
  int tnv = m->getDownward(t, 0, tvs);
  PCU_ALWAYS_ASSERT(tnv == 4 && fnv == 3);
  for (int i = 0; i < tnv; i++) {
    if (apf::findIn(fvs, fnv, tvs[i]) == -1)
      return tvs[i];
  }
  return 0;
}

/*
 * Given a face and one of its edges, the vertex of the face
 * opposite to the given edge is returned.
 */
static apf::MeshEntity* getFaceOppVert(
    apf::Mesh* m, apf::MeshEntity* f, apf::MeshEntity* e)
{
  apf::Downward evs;
  int env = m->getDownward(e, 0, evs);
  apf::Downward fvs;
  int fnv = m->getDownward(f, 0, fvs);
  PCU_ALWAYS_ASSERT(env == 2 && fnv == 3);
  for (int i = 0; i < fnv; i++) {
    if (apf::findIn(evs, env, fvs[i]) == -1)
      return fvs[i];
  }
  return 0;
}

static apf::Vector3 computeFaceNormal(
    apf::Mesh* m, apf::MeshEntity* f, apf::Vector3 const& p)
{
  // Compute face normal using face Jacobian,
  // so it can also be used for curved mesh elements
  apf::MeshElement* me = apf::createMeshElement(m, f);
  apf::Matrix3x3 J;
  apf::getJacobian(me, p, J);
  apf::destroyMeshElement(me);

  apf::Vector3 g1 = J[0];
  apf::Vector3 g2 = J[1];
  apf::Vector3 n = apf::cross( g1, g2 );
  return n.normalize();
}

apf::Vector3 computeFaceOutwardNormal(apf::Mesh* m,
    apf::MeshEntity* t, apf::MeshEntity* f, apf::Vector3 const& p)
{
  apf::Vector3 n = computeFaceNormal(m, f, p);

  // orient the normal outwards from the tet
  apf::MeshEntity* oppVert = getTetOppVert(m, t, f);
  apf::Vector3 vxi = apf::Vector3(0.,0.,0.);

  apf::Vector3 txi;
  m->getPoint(oppVert, 0, txi);

  apf::MeshElement* fme = apf::createMeshElement(m, f);
  apf::Vector3 pxi;
  apf::mapLocalToGlobal(fme, p, pxi);
  apf::destroyMeshElement(fme);

  apf::Vector3 pxiTotxi = txi - pxi;
  if (pxiTotxi*n > 0) {
    n = n*-1.;
  }
  return n;
}

/*
 * computes local flux integral restricted to
 * an edge of a tet element
 */
static double getLocalFluxIntegral(EdgePatch* ep, apf::MeshEntity* tet)
{
  PCU_ALWAYS_ASSERT(ep->mesh->getType(tet) == apf::Mesh::TET);

  double fluxIntegral = 0.0;
  // 1. get faces of the tet in the patch
  apf::Downward f;
  int nf = ep->mesh->getDownward(tet, 2, f);
  std::vector<apf::MeshEntity*> patchFaces;
  for (int i = 0; i < nf; i++) {
    if(std::find(ep->faces.begin(), ep->faces.end(), f[i]) != ep->faces.end())
      patchFaces.push_back(f[i]);
  }
  PCU_ALWAYS_ASSERT(patchFaces.size() == 2);

  //  2. loop over the patch faces
  for (unsigned int i = 0; i < patchFaces.size(); i++) {
    double fluxFaceIntegral = 0.0;
    // 3. get upward tets of the current face
    apf::Up up;
    apf::MeshEntity* currentFace = patchFaces[i];
    ep->mesh->getUp(currentFace, up);
    if (crv::isBoundaryEntity(ep->mesh, currentFace))
      PCU_ALWAYS_ASSERT( up.n == 1);
    else
      PCU_ALWAYS_ASSERT( up.n == 2);

    apf::MeshEntity* firstTet = up.e[0];
    apf::MeshEntity* secondTet;
    if (up.n == 2)
      secondTet  = up.e[1];

    // 4. findIn edge in downward edges of current face
    apf::Downward e;
    int ne = ep->mesh->getDownward(currentFace, 1, e);
    PCU_ALWAYS_ASSERT(ne == 3);
    int ei = apf::findIn(e, ne, ep->entity);

    // 5. count integration points for flux face integral
    apf::FieldShape* fs = ep->equilibration->ef->getShape();
    int int_order = 2 * fs->getOrder();
    apf::MeshElement* fme = apf::createMeshElement(ep->mesh, currentFace);
    apf::Element* fel = apf::createElement(ep->equilibration->ef, fme);
    int np = apf::countIntPoints(fme, int_order);

    // loop over integration points
    apf::Vector3 p, tet1xi, tet2xi, curl1, curl2, curl,
      fn1, fn2, tk, vshape;
    for (int n = 0; n < np; n++) {
      apf::getIntPoint(fme, int_order, n, p);
      double weight = apf::getIntWeight(fme, int_order, n);
      apf::Matrix3x3 fJ;
      apf::getJacobian(fme, p, fJ);
      double jdet = apf::getJacobianDeterminant(
          fJ, apf::getDimension(ep->mesh, currentFace));

      // compute face outward normals wrt tets
      if (tet == firstTet)
        fn1 = computeFaceOutwardNormal(ep->mesh, firstTet, currentFace, p);
      else
        fn1 = computeFaceOutwardNormal(ep->mesh, secondTet, currentFace, p);
      if (up.n == 2) {
        if (tet == firstTet)
          fn2 = computeFaceOutwardNormal(ep->mesh, secondTet, currentFace, p);
        else
          fn2 = computeFaceOutwardNormal(ep->mesh, firstTet, currentFace, p);
      }

      curl.zero();
      // compute curl1
      tet1xi = apf::boundaryToElementXi(ep->mesh, currentFace, firstTet, p);
      apf::MeshElement* me1 = apf::createMeshElement(ep->mesh, firstTet);
      apf::Element* el1 = apf::createElement(ep->equilibration->ef, me1);
      apf::getCurl(el1, tet1xi, curl1);
      apf::Vector3 temp1 = apf::cross(fn1, curl1);
      curl += temp1;
      apf::destroyElement(el1);
      apf::destroyMeshElement(me1);

      // compute curl2
      if (up.n == 2) {
        tet2xi = apf::boundaryToElementXi(ep->mesh, currentFace, secondTet, p);
        apf::MeshElement* me2 = apf::createMeshElement(ep->mesh, secondTet);
        apf::Element* el2 = apf::createElement(ep->equilibration->ef, me2);
        apf::getCurl(el2, tet2xi, curl2);
        apf::Vector3 temp2 = apf::cross(fn2, curl2);
        curl += (temp2 * -1.);
        curl = curl * 1./2.;
        apf::destroyElement(el2);
        apf::destroyMeshElement(me2);
      }

      // compute tk (inter-element averaged flux)
      tk = curl;

      // compute vector shape
      int type = apf::Mesh::TRIANGLE;
      int nd = apf::countElementNodes(fs, type);
      apf::NewArray<apf::Vector3> vectorshapes (nd);
      apf::getVectorShapeValues(fel, p, vectorshapes);
      vshape = vectorshapes[ei];

      // TODO why?
      int which, rotate; bool flip;
      apf::getAlignment(
          ep->mesh, currentFace, ep->entity, which, flip, rotate);
      if (flip) {
      vshape = vshape * -1.;
      }

      // compute integral
      fluxFaceIntegral += (tk * vshape) * weight * jdet;
    }
    apf::destroyElement(fel);
    apf::destroyMeshElement(fme);

    fluxIntegral += fluxFaceIntegral;
  }
  return fluxIntegral;
}

static void assembleEdgePatchRHS(EdgePatch* p)
{
  if (p->isOnBdry) {
    p->b.resize(p->tets.size() + p->faces.size());
    p->b.zero();
  }
  else {
    p->b.resize(p->tets.size());
    p->b.zero();
  }
  int ne = p->tets.size();
  int nf = p->faces.size();
  for (int i = 0; i < ne; i++) {
    apf::MeshEntity* tet = p->tets[i];
    double blfIntegral = getLocalEdgeBLF(p, tet);
    double lfIntegral = getLocalEdgeLF(p, tet);
    double fluxIntegral = getLocalFluxIntegral(p, tet);
    if(p->isOnBdry)
      p->b(nf+i) = blfIntegral - lfIntegral - fluxIntegral;
    else
      p->b(i) = blfIntegral - lfIntegral - fluxIntegral;
  }
}

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

/*
 * Orders tets and faces in an edge cavity in a
 * clockwise (or ccw) direction around the edge.
 */
static void getOrderedTetsandFaces(apf::Mesh* mesh, apf::MeshEntity* edge,
     EntityVector& tets, EntityVector& faces)
{
  tets.clear();
  faces.clear();
  if( ! crv::isBoundaryEntity(mesh, edge) ) {
    apf::MeshEntity* currentFace = mesh->getUpward(edge, 0);
    apf::Up up;
    mesh->getUp(currentFace, up);
    PCU_ALWAYS_ASSERT(up.n == 2);
    apf::MeshEntity* firstTet = up.e[0];
    apf::MeshEntity* nextTet  = up.e[1];
    tets.push_back(firstTet);
    apf::MeshEntity* firstFace = getTetOppFaceSharingEdge(mesh, firstTet,
        currentFace, edge);
    faces.push_back(firstFace);

    while (nextTet != firstTet) {
      tets.push_back(nextTet);
      faces.push_back(currentFace);
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
  else {
    apf::Up up;
    mesh->getUp(edge, up);
    apf::MeshEntity* firstFace;
    for (int i = 0; i < up.n; i++) {
      if ( crv::isBoundaryEntity(mesh, up.e[i]) ) {
        firstFace = up.e[i]; break;
      }
    }
    faces.push_back(firstFace);
    mesh->getUp(firstFace, up);
    PCU_ALWAYS_ASSERT(up.n == 1);
    apf::MeshEntity* firstTet = up.e[0];
    tets.push_back(firstTet);

    apf::MeshEntity* nextFace = getTetOppFaceSharingEdge(mesh, firstTet,
        firstFace, edge);
    apf::MeshEntity* nextTet = firstTet;
    mesh->getUp(nextFace, up);
    while( up.n == 2) {
      faces.push_back(nextFace);
      if (nextTet != up.e[0])
        nextTet = up.e[0];
      else
        nextTet = up.e[1];
      tets.push_back(nextTet);

      nextFace = getTetOppFaceSharingEdge(mesh, nextTet, nextFace, edge);
      mesh->getUp(nextFace, up);
    }
    faces.push_back(nextFace);
  }
}

void getClockwiseTetsandFaces(EdgePatch* p)
{
  if (p->isOnBdry) {
    if (p->tets.size() > 1) {
      apf::MeshEntity* secondFace = p->faces[1];
      apf::MeshEntity* firstTet = p->tets[0];
      apf::MeshEntity* secondTet = p->tets[1];

      apf::Downward firstTetFaces, secondTetFaces;
      int n_firstTetFaces = p->mesh->getDownward(firstTet, 2, firstTetFaces);
      int n_secondTetFaces = p->mesh->getDownward(secondTet, 2, secondTetFaces);
      int fi1 = apf::findIn(firstTetFaces, n_firstTetFaces, secondFace);
      int fi2 = apf::findIn(secondTetFaces, n_secondTetFaces, secondFace);
      PCU_ALWAYS_ASSERT(fi1 != -1 && fi2 != -1);

      // first tet opp vertex crd
      apf::MeshEntity* firstTetOppVert = getTetOppVert(
        p->mesh, firstTet, secondFace);
      apf::Vector3 firstTetOppVertCrd;
      p->mesh->getPoint(firstTetOppVert, 0, firstTetOppVertCrd);

      // second tet opp vertex crd
      apf::MeshEntity* secondTetOppVert = getTetOppVert(
        p->mesh, secondTet, secondFace);
      apf::Vector3 secondTetOppVertCrd;
      p->mesh->getPoint(secondTetOppVert, 0, secondTetOppVertCrd);

      // normal to the face
      apf::Downward edge_vertices;
      p->mesh->getDownward(p->entity, 0, edge_vertices);
      apf::Vector3 p0, p1, p2;
      p->mesh->getPoint(edge_vertices[0], 0, p0);
      p->mesh->getPoint(edge_vertices[1], 0, p1);
      apf::MeshEntity* faceOppVert = getFaceOppVert(
        p->mesh, secondFace, p->entity);
      p->mesh->getPoint(faceOppVert, 0, p2);

      apf::Vector3 normal = apf::cross(p1-p0, p2-p0);

      // direction vectors from p0 to opp tet verts
      apf::Vector3 vFirst = firstTetOppVertCrd - p0;
      apf::Vector3 vLast = secondTetOppVertCrd - p0;

      if ((vFirst * normal > 0)  && (vLast * normal < 0)) {
        // reverse list of tets and faces
        std::reverse(p->tets.begin(), p->tets.end());
        std::reverse(p->faces.begin(), p->faces.begin());
      }
      else if ((vFirst * normal < 0)  && (vLast * normal > 0)) {
        cout << "already clockwise" << endl;
      }
      else if ((vFirst * normal < 0)  && (vLast * normal < 0)) {
        cout << "failed clockwise. ABORT" << endl;
      }
      else if ((vFirst * normal > 0)  && (vLast * normal > 0)) {
        cout << "failed clockwise. ABORT" << endl;
      }
    }
  }
  else {
    apf::MeshEntity* firstFace = p->faces[0];
    apf::MeshEntity* firstTet = p->tets[0];
    apf::MeshEntity* lastTet = p->tets[p->tets.size()-1];

    apf::Downward firstTetFaces, lastTetFaces;
    int n_firstTetFaces = p->mesh->getDownward(firstTet, 2, firstTetFaces);
    int n_lastTetFaces = p->mesh->getDownward(lastTet, 2, lastTetFaces);
    int fi1 = apf::findIn(firstTetFaces, n_firstTetFaces, firstFace);
    int filast = apf::findIn(lastTetFaces, n_lastTetFaces, firstFace);
    PCU_ALWAYS_ASSERT(fi1 != -1 && filast != -1);

    // first tet opp vertex crd
    apf::MeshEntity* firstTetOppVert = getTetOppVert(
        p->mesh, firstTet, firstFace);
    apf::Vector3 firstTetOppVertCrd;
    p->mesh->getPoint(firstTetOppVert, 0, firstTetOppVertCrd);

    // last tet opp vertex crd
    apf::MeshEntity* lastTetOppVert = getTetOppVert(
        p->mesh, lastTet, firstFace);
    apf::Vector3 lastTetOppVertCrd;
    p->mesh->getPoint(lastTetOppVert, 0, lastTetOppVertCrd);

    // normal to the face
    apf::Downward edge_vertices;
    p->mesh->getDownward(p->entity, 0, edge_vertices);
    apf::Vector3 p0, p1, p2;
    p->mesh->getPoint(edge_vertices[0], 0, p0);
    p->mesh->getPoint(edge_vertices[1], 0, p1);
    apf::MeshEntity* faceOppVert = getFaceOppVert(
        p->mesh, firstFace, p->entity);
    p->mesh->getPoint(faceOppVert, 0, p2);

    apf::Vector3 normal = apf::cross(p1-p0, p2-p0);

    // direction vectors from p0 to opp tet verts
    apf::Vector3 vFirst = firstTetOppVertCrd - p0;
    apf::Vector3 vLast = lastTetOppVertCrd - p0;

    if ((vFirst * normal > 0)  && (vLast * normal < 0)) {
      // reverse list of tets and faces
      std::reverse(p->tets.begin(), p->tets.end());
      std::reverse(p->faces.begin(), p->faces.begin());
    }
    else if ((vFirst * normal < 0)  && (vLast * normal > 0)) {
    }
    else if ((vFirst * normal < 0)  && (vLast * normal < 0)) {
      cout << "failed clockwise. ABORT" << endl;
    }
    else if ((vFirst * normal > 0)  && (vLast * normal > 0)) {
      cout << "failed clockwise. ABORT" << endl;
    }
  }
}

static void runErm(EdgePatch* ep)
{
  getOrderedTetsandFaces(ep->mesh, ep->entity, ep->tets, ep->faces);
  //getClockwiseTetsandFaces(ep);
  assembleEdgePatchLHS(ep);
  assembleEdgePatchRHS(ep);
  mth::solveFromQR(ep->qr.Q, ep->qr.R, ep->b, ep->x);

  if (!ep->isOnBdry) { // solve At*mu = g for g
    mth::Vector<double> temp(ep->tets.size());
    mth::multiply(ep->At, ep->x, temp);
    for (size_t i = 0; i < ep->tets.size(); i++) {
      ep->x(i) = temp(i);
    }
  }

  int nf = ep->faces.size();
  for(int i = 0; i < nf; i++) {
    apf::MeshEntity* face = ep->faces[i];

    apf::Downward e;
    int ned = ep->mesh->getDownward(face, 1, e);
    int ei = apf::findIn(e, ned, ep->entity);
    PCU_ALWAYS_ASSERT(ned == 3 && ei != -1);

    double components[3];
    apf::getComponents(ep->equilibration->g, face, 0, components);
    components[ei] = ep->x(i);
    apf::setComponents(ep->equilibration->g, face, 0, components);
  }
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

apf::Field* equilibrateResiduals(apf::Field* f)
{
  apf::Field* g = createPackedField(
      apf::getMesh(f), "g", 3, apf::getConstant(2) );
  Equilibration equilibration;
  setupEquilibration(&equilibration, f, g);
  EdgePatchOp op(&equilibration);
  op.applyToDimension(1); // edges
  return g;
}


}

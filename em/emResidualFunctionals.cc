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
  /* input scalar field containing Nedelec dofs for electric field */
  apf::Field* ef;
  /* output field containing correction values */
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
  apf::MeshIterator* it = apf::getMesh(f)->begin(2);
  while ((face = apf::getMesh(f)->iterate(it))) {
    apf::setComponents(eq->g, face, 0, zeros);
  }
  apf::getMesh(f)->end(it);
}

struct QRDecomp {
  mth::Matrix<double> Q;
  mth::Matrix<double> R;
};

typedef std::vector<apf::MeshEntity*> EntityVector;

struct EdgePatch {
  apf::Mesh* mesh;
  Equilibration* equilibration;
  /* the entity around which the patch
     is centered. a patch collects entities
     (faces & tets) around this edge entity */
  apf::MeshEntity* entity;
  bool isOnBdry;
  EntityVector tets;
  EntityVector faces;
  mth::Matrix<double> A;
  mth::Vector<double> b;
  mth::Vector<double> x;
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
    p->faces.push_back(e);
  if(p->mesh->getType(e) == apf::Mesh::TET)
    p->tets.push_back(e);
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
  // REMOVE ALL
  int ne = adjacent.size();
  std::cout << " Unordered Tets, ne: " << ne << std::endl; // REMOVE
  std::cout << " Center of Unordered tets" << std::endl;
  for (int i = 0; i < ne; i++) {
    apf::MeshEntity* tet = adjacent[i];;
    apf::Vector3 center = apf::getLinearCentroid(p->mesh, tet);
    std::cout << i << ":  " << center << std::endl;
  } //////////



  p->mesh->getAdjacent(p->entity, 2, adjacent);
  addEntitiesToPatch(p, adjacent);

  // REMOVE ALL
  int nf = adjacent.size();
  std::cout << " Unordered Faces, nf: " << nf << std::endl; // REMOVE
  std::cout << " Center of Unordered Faces" << std::endl;
  for (int i = 0; i < nf; i++) {
    apf::MeshEntity* face = adjacent[i];
    apf::Vector3 center = apf::getLinearCentroid(p->mesh, face);
    std::cout << i << ":  " << center << std::endl;
  }
  std::cout << "----------------" << std::endl;

  return true;
}

static bool buildEdgePatch(EdgePatch* p, apf::CavityOp* o)
{
  if (!getInitialEdgePatch(p, o)) return false;
  return true;
}

static void assembleEdgePatchLHS(EdgePatch* p)
{
  int ne = p->tets.size();
  int nf = p->faces.size();
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
  }
  mth::decomposeQR(p->A, p->qr.Q, p->qr.R);
}

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
  apf::destroyElement(el);
  apf::destroyMeshElement(me);
}

void assembleVectorMassElementMatrix(apf::Mesh* mesh, apf::MeshEntity* e,
    apf::Field* f, mth::Matrix<double>& elmat)
{
  apf::FieldShape* fs = f->getShape();
  int type = mesh->getType(e);
  PCU_ALWAYS_ASSERT(type == apf::Mesh::TET || type == apf::Mesh::TRIANGLE);
  int nd = apf::countElementNodes(fs, type);
  int dim = apf::getDimension(mesh, e);
  double w;

  apf::NewArray<apf::Vector3> vectorshapes(nd);
  elmat.resize(nd,nd);

  apf::MeshElement* me = apf::createMeshElement(mesh, e);
  apf::Element* el = apf::createElement(f, me);
  int int_order = 2 * fs->getOrder();
  int np = apf::countIntPoints(me, int_order); // int points required

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
    mth::Matrix<double> vectorShapes (nd, dim); // TODO opt clean
    for (int j = 0; j < nd; j++)
      for (int k = 0; k < dim; k++)
        vectorShapes(j,k) = vectorshapes[j][k];

    mth::Matrix<double> vectorShapesT (dim, nd);
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

// computes local bilinear form integral restricted to
// an edge of a tet element
static double getLocalEdgeBLF(EdgePatch* ep, apf::MeshEntity* tet)
{
  /* findIn edge in downward edges of tet
  apf::Downward e;
  int ne = p->mesh->getDownward(tet, 1, e);
  PCU_ALWAYS_ASSERT(ne == 6);
  int ei = apf::findIn(e, ne, p->entity);
  // get Element Dofs
  apf::MeshElement* me = apf::createMeshElement(p->mesh, tet);
  apf::Element* el = apf::createElement(p->equilibration->ef, me);
  int type = p->mesh->getType(tet);
  int nd = apf::countElementNodes(el->getFieldShape(), type);
  apf::NewArray<double> d (nd);
  el-getElementDofs(d);
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
  mth::Vector<double> blf_integrals (nd);
  mth::multiply(elmat, dofs, blf_integrals);

  apf::destroyElement(el);
  apf::destroyMeshElement(me);*/

  // pick edge index from the resulting vector
  // negation of negative ND dofs
  /*int which, rotate; bool flip;
  apf::getAlignment(p->mesh, tet, p->entity, which, flip, rotate);
  if (flip)
    blf_integrals(ei) = -1*blf_integrals(ei);*/

  /*cout << "local blf " << blf_integrals(ei) << endl;
  return blf_integrals(ei);*/
  // 0. findIn edge in downward edges of tet
  apf::Downward e;
  int ne = ep->mesh->getDownward(tet, 1, e);
  PCU_ALWAYS_ASSERT(ne == 6);
  int ei = apf::findIn(e, ne, ep->entity);


  apf::FieldShape* fs = ep->equilibration->ef->getShape();
  int type = ep->mesh->getType(tet);
  int dim = apf::getDimension(ep->mesh, tet);


  apf::MeshElement* me = apf::createMeshElement(ep->mesh, tet);
  apf::Element* el = apf::createElement(ep->equilibration->ef, me);
  int nd = apf::countElementNodes(el->getFieldShape(), type);
  int int_order = 2 * fs->getOrder();
  int np = apf::countIntPoints(me, int_order);
  double w;


  apf::NewArray<apf::Vector3> curlshapes(nd);
  apf::NewArray<apf::Vector3> vectorshapes(nd);
  mth::Matrix<double> phys_curlshapes(nd, dim);  // TODO clean


  // 1. Curl Curl Integration
  double curlcurl = 0.0;

  apf::Vector3 p, curl;
  for (int i = 0; i < np; i++) {
    apf::getIntPoint(me, int_order, i, p);
    double weight = apf::getIntWeight(me, int_order, i);
    apf::Matrix3x3 J;
    apf::getJacobian(me, p, J);
    //double jdet = apf::getJacobianDeterminant(J, dim);
    w = weight; // TODO check why do not need division by jdet

    // get curl vector
    apf::getCurl(el, p, curl);

    // get curlshape values // TODO CLEAN use getCurlShapeValues
    el->getShape()->getLocalVectorCurls(ep->mesh, tet, p, curlshapes);
    phys_curlshapes.zero();
    for (int i = 0; i < nd; i++)
      for (int j = 0; j < dim; j++)
        for (int k = 0; k < dim; k++)
          phys_curlshapes(i,j) += curlshapes[i][k] * J[k][j];

    // pick the ei curlshape // TODO Clean apf::Vector
    apf::Vector3 cshape;
    cshape[0] = phys_curlshapes(ei, 0);
    cshape[1] = phys_curlshapes(ei, 1);
    cshape[2] = phys_curlshapes(ei, 2);

    // multiply
    curlcurl += (curl * cshape) * w;
  }

  // 2. Vector Mass Integration
  int_order = 2 * fs->getOrder();
  np = apf::countIntPoints(me, int_order); // int points required

  double vectormass = 0.0;
  apf::Vector3 vvalue;
  for (int i = 0; i < np; i++) {
    apf::getIntPoint(me, int_order, i, p);
    double weight = apf::getIntWeight(me, int_order, i);
    apf::Matrix3x3 J;
    apf::getJacobian(me, p, J);
    double jdet = apf::getJacobianDeterminant(J, dim);

    w = weight * jdet;

    apf::getVector(el, p, vvalue);

    apf::getVectorShapeValues(el, p, vectorshapes);
    apf::Vector3 vshape = vectorshapes[ei];

    vectormass += (vvalue * vshape) * w;
  }

  double blf_integral = curlcurl + vectormass;

  int which, rotate; bool flip;
  apf::getAlignment(ep->mesh, tet, ep->entity, which, flip, rotate);
  if (flip)
    blf_integral = -1*blf_integral;
  cout << "local blf " << blf_integral << endl;
  return blf_integral;
}


// TODO QUESTION redo this to allow user access from outside
void pumiUserFunction(const apf::Vector3& x, mth::Vector<double>& f,
    apf::MeshEntity* e, apf::Mesh* mesh)
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
    pumiUserFunction(global, val, e, apf::getMesh(f)); // eval vector function
    val *= w;

    mth::Matrix<double> vectorShapes (nd, dim);
    for (int i = 0; i < nd; i++)
      for (int j = 0; j < dim; j++)
        vectorShapes(i,j) = vectorshapes[i][j];
    mth::Vector<double> V (nd);
    V.zero();
    mth::multiply(vectorShapes, val, V);
    elvect += V;
  }
  apf::destroyElement(el);
  apf::destroyMeshElement(me);
}

// computes local linear form integral restricted to
// an edge of a tet element
static double getLocalEdgeLF(EdgePatch* p, apf::MeshEntity* tet)
{
  // findIn edge in downward edges of tet
  apf::Downward e;
  int ne = p->mesh->getDownward(tet, 1, e);
  PCU_ALWAYS_ASSERT(ne == 6);
  int ei = apf::findIn(e, ne, p->entity);
  // assemble Domain LF Vector
  mth::Vector<double> elvect;
  assembleDomainLFElementVector(p->mesh, tet,
      p->equilibration->ef, elvect);
  int which, rotate; bool flip;
  apf::getAlignment(p->mesh, tet, p->entity, which, flip, rotate);
  if (flip)
    elvect(ei) = -1*elvect(ei);
  cout << "local lf " << elvect(ei) << endl;
  return elvect(ei);
}

// Given a tet and one of its faces, the vertex of the tet
// opposite to the face is returned.
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

static apf::Vector3 computeFaceNormal(apf::Mesh* m,
    apf::MeshEntity* f, apf::Vector3 const& p)
{
  // Compute face normal using face Jacobian
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

  // Orient the normal outwards from the tet
  apf::MeshEntity* oppVert = getTetOppVert(m, t, f);
  apf::Vector3 vxi = apf::Vector3(0.,0.,0.);

  // get global coordinates of the vertex
  apf::Vector3 txi;
  m->getPoint(oppVert, 0, txi);

  // get global coordinates of the point on the face
  apf::MeshElement* fme = apf::createMeshElement(m, f);
  apf::Vector3 pxi;
  apf::mapLocalToGlobal(fme, p, pxi);
  apf::destroyMeshElement(fme);

  apf::Vector3 pxiTotxi = txi - pxi;
  //std::cout << "dot product " << pxiTotxi*n << std::endl;
  if (pxiTotxi*n > 0) {
    n = n*-1.;
  }
  return n;
}

// computes local flux integral restricted to
// an edge of a tet element
static double getLocalFluxIntegral(EdgePatch* ep, apf::MeshEntity* tet)
{
  double fluxIntegral = 0.0;
  // 1. get faces of the tet in the patch
  apf::Downward f;
  int nf = ep->mesh->getDownward(tet, 2, f);
  PCU_ALWAYS_ASSERT(nf == 4);
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
    std::cout << "np " << np << std::endl; // REMOVE

    // loop over integration points
    apf::Vector3 p, tet1xi, tet2xi, curl1, curl2, curl,
      fnormal1, fnormal2, tk, vshape;
    for (int n = 0; n < np; n++) {
      apf::getIntPoint(fme, int_order, n, p);
      double weight = apf::getIntWeight(fme, int_order, n);
      apf::Matrix3x3 fJ;
      apf::getJacobian(fme, p, fJ);
      double jdet = apf::getJacobianDeterminant(
          fJ, apf::getDimension(ep->mesh, currentFace));

      // compute face outward normals wrt tets
      if (tet == firstTet)
        fnormal1 = computeFaceOutwardNormal(ep->mesh, firstTet, currentFace, p);
      else
        fnormal1 = computeFaceOutwardNormal(ep->mesh, secondTet, currentFace, p);
      if (up.n == 2) {
        if (tet == firstTet)
          fnormal2 = computeFaceOutwardNormal(ep->mesh, secondTet, currentFace, p);
        else
          fnormal2 = computeFaceOutwardNormal(ep->mesh, firstTet, currentFace, p);
        std::cout << "normal1 " << fnormal1 << std::endl; // REMOVE
        std::cout << "normal2 " << fnormal2 << std::endl; // REMOVE
      }

      curl.zero();
      // compute curl1
      tet1xi = apf::boundaryToElementXi(ep->mesh, currentFace, firstTet, p);
      apf::MeshElement* me1 = apf::createMeshElement(ep->mesh, firstTet);
      apf::Element* el1 = apf::createElement(ep->equilibration->ef, me1);
      apf::getCurl(el1, tet1xi, curl1);
      apf::Vector3 temp1 = apf::cross(fnormal1, curl1);
      curl += temp1;
      apf::destroyElement(el1);
      apf::destroyMeshElement(me1);

      // compute curl2
      if (up.n == 2) {
        tet2xi = apf::boundaryToElementXi(ep->mesh, currentFace, secondTet, p);
        apf::MeshElement* me2 = apf::createMeshElement(ep->mesh, secondTet);
        apf::Element* el2 = apf::createElement(ep->equilibration->ef, me2);
        apf::getCurl(el2, tet2xi, curl2);
        apf::Vector3 temp2 = apf::cross(fnormal2, curl2);
        curl += (temp2 * -1.);
        apf::destroyElement(el2);
        apf::destroyMeshElement(me2);
      }

      // compute tk (inter-element averaged flux)
      //tk = apf::cross(fnormal, curl);
      tk = curl * 1./2.;
      std::cout << "tk " << tk << std::endl; // REMOVE

      // compute vector shape
      int type = apf::Mesh::TRIANGLE;
      int nd = apf::countElementNodes(fs, type);
      apf::NewArray<apf::Vector3> vectorshapes (nd);
      apf::getVectorShapeValues(fel, p, vectorshapes);
      vshape = vectorshapes[ei];

      // compute integral
      fluxFaceIntegral += (tk * vshape) * weight * jdet;
    }
    apf::destroyElement(fel);
    apf::destroyMeshElement(fme);
    fluxIntegral += fluxFaceIntegral;
    std::cout << "flux Face integral " << fluxFaceIntegral << std::endl; // REMOVE
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
  double testblflf = 0.0; // REMOVE
  double testflux = 0.0; // REMOVE
  int ne = p->tets.size();
  int nf = p->faces.size();
  for (int i = 0; i < ne; i++) {
    apf::MeshEntity* tet = p->tets[i];
    double blfIntegral = getLocalEdgeBLF(p, tet);
    double lfIntegral = getLocalEdgeLF(p, tet);
    double fluxIntegral = getLocalFluxIntegral(p, tet);
    testblflf += (blfIntegral - lfIntegral);
    testflux += fluxIntegral;
    cout << "Blf Integral" << blfIntegral << endl;
    cout << "Lf Integral" << lfIntegral << endl;
    if(p->isOnBdry)
      p->b(nf+i) = blfIntegral - lfIntegral - fluxIntegral;
    else
      p->b(i) = blfIntegral - lfIntegral - fluxIntegral;
  }
  std::cout << "Blf - lf integrals = " << testblflf << std::endl;
  std::cout << "Flux Integral Sum  = " << testflux << std::endl;
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
     EntityVector& tets, EntityVector& faces)
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
  else { // boundary edge patches
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

static void runErm(EdgePatch* p)
{
  getOrderedTetsandFaces(p->mesh, p->entity, p->tets, p->faces);
  assembleEdgePatchLHS(p);
  assembleEdgePatchRHS(p);
  mth::solveFromQR(p->qr.Q, p->qr.R, p->b, p->x);

  bool debug = true; // REMOVE DEBUG
  if (debug) {
    if (p->isOnBdry)
      cout << "Boundary Patch" << endl;
    else
      cout << "Interior Patch" << endl;
    cout << "ne " << p->tets.size() << " nf " << p->faces.size() << endl;
    cout << "LHS Matrix" << endl;
    cout << p->A << endl;

    cout << "RHS Vector" << endl;
    cout << p->b << endl;
    double rhs_sum = 0.0;
    for(unsigned int i = 0; i < p->b.size(); i++)
      rhs_sum += p->b(i);
    cout << "rhs_sum " << rhs_sum << endl;
    if ( (abs(rhs_sum) > 1e-8) && (!p->isOnBdry) )
      cout << "nonzero! error" << endl;
  }

  int nf = p->faces.size();
  for(int i = 0; i < nf; i++) {
    // get face entitiy
    apf::MeshEntity* face = p->faces[i];
    // get ei of edge in face downward edges
    apf::Downward e;
    int ned = p->mesh->getDownward(face, 1, e);
    int ei = apf::findIn(e, ned, p->entity);
    PCU_ALWAYS_ASSERT(ned == 3 && ei != -1);
    // save the g value at that index on that face in the g field
    double components[3];
    apf::getComponents(p->equilibration->g, face, 0, components);
    components[ei] = p->x(i);
    apf::setComponents(p->equilibration->g, face, 0, components);
  }

  // REMOVE ALL BELOW
  int ne = p->tets.size();
  //int nf = p->faces.size();
  std::cout << " Reordered Tets, ne: " << ne << " nf " << nf << std::endl; // REMOVE
  std::cout << " Center of Reordered tets" << std::endl;
  for (int i = 0; i < ne; i++) {
    apf::MeshEntity* tet = p->tets[i];
    apf::Vector3 center = apf::getLinearCentroid(p->mesh, tet);
    std::cout << i << ":  " << center << std::endl;
  }
  std::cout << " Center of Reordered Faces" << std::endl;
  for (int i = 0; i < nf; i++) {
    apf::MeshEntity* face = p->faces[i];
    apf::Vector3 center = apf::getLinearCentroid(p->mesh, face);
    std::cout << i << ":  " << center << std::endl;
  }
  std::cout << "----------------" << std::endl;

  std::cout << "x" << std::endl; // REMOVE
  std::cout << p->x << std::endl; // REMOVE
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
  apf::Field* g = createPackedField( apf::getMesh(f), "g", 3, apf::getConstant(2) );
  Equilibration equilibration;
  setupEquilibration(&equilibration, f, g);
  EdgePatchOp op(&equilibration);
  op.applyToDimension(1); // edges
  return g;
}







}

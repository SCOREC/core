/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include <iostream>
#include <cstdlib>

#include "apfElement.h"
#include "crv.h"
#include "crvShape.h"

#include "em.h"
using namespace std;
namespace em {

static void computeResidualBLF(apf::Mesh* mesh, apf::MeshEntity* e,
  apf::Field* f, apf::Field* fp1, mth::Vector<double>& blf)
{
  apf::FieldShape* fp1s = fp1->getShape();
  int type = mesh->getType(e);
  PCU_ALWAYS_ASSERT(type == apf::Mesh::TET);
  int nd = apf::countElementNodes(fp1s, type);
  int dim = apf::getDimension(mesh, e);
  double w;

  apf::NewArray<apf::Vector3> curlshape(nd);
  apf::NewArray<apf::Vector3> vectorshape(nd);
  mth::Matrix<double> phys_curlshape(nd, dim);
  mth::Matrix<double> vectorShape(nd, dim);
  mth::Vector<double> curlcurl_vec (nd);
  mth::Vector<double> mass_vec (nd);
  blf.resize(nd);

  apf::MeshElement* me = apf::createMeshElement(mesh, e);
  apf::Element* fp1el = apf::createElement(fp1, me);
  apf::Element* fel = apf::createElement(f, me);
  int int_order = 2 * fp1s->getOrder();
  int np = apf::countIntPoints(me, int_order);

  // 1. Compute Curl Curl Integration
  curlcurl_vec.zero();
  apf::Vector3 p, curl;
  for (int i = 0; i < np; i++) {
    apf::getIntPoint(me, int_order, i, p);
    double weight = apf::getIntWeight(me, int_order, i);
    apf::Matrix3x3 J;
    apf::getJacobian(me, p, J);
    //double jdet = apf::getJacobianDeterminant(J, dim);
    w = weight; // TODO check why do not need division by jdet

    // get curl vector
    apf::getCurl(fel, p, curl);

    // get curlshape values // TODO CLEAN use getCurlShapeValues
    fp1el->getShape()->getLocalVectorCurls(mesh, e, p, curlshape);
    phys_curlshape.zero();
    for (int i = 0; i < nd; i++)
      for (int j = 0; j < dim; j++)
        for (int k = 0; k < dim; k++)
          phys_curlshape(i,j) += curlshape[i][k] * J[k][j];

    // multiply
    mth::Vector<double> V (nd);
    V.zero();
    mth::Vector<double> c (dim);
    c(0) = curl[0]; c(1) = curl[1]; c(2) = curl[2];
    mth::multiply(phys_curlshape, c, V);
    V *= w;

    curlcurl_vec += V;
  }

  // 2. Compute Vector Mass Integration
  int_order = 2 * fp1s->getOrder();
  np = apf::countIntPoints(me, int_order); // int points required

  mass_vec.zero();
  apf::Vector3 vvalue;
  for (int i = 0; i < np; i++) {
    apf::getIntPoint(me, int_order, i, p);
    double weight = apf::getIntWeight(me, int_order, i);
    apf::Matrix3x3 J;
    apf::getJacobian(me, p, J);
    double jdet = apf::getJacobianDeterminant(J, dim);

    w = weight * jdet;

    apf::getVector(fel, p, vvalue);

    apf::getVectorShapeValues(fp1el, p, vectorshape);
    vectorShape.zero();
    for (int i = 0; i < nd; i++)
      for (int j = 0; j < dim; j++)
        vectorShape(i,j) = vectorshape[i][j];


    mth::Vector<double> V (nd);
    V.zero();
    mth::Vector<double> v (dim);
    v(0) = vvalue[0]; v(1) = vvalue[1]; v(2) = vvalue[2];
    mth::multiply(vectorShape, v, V);
    V *= w;

    mass_vec += V;
  }

  // 3. get result
  blf.zero();
  blf += curlcurl_vec;
  blf += mass_vec;

  apf::destroyMeshElement(me);
  apf::destroyElement(fp1el);
  apf::destroyElement(fel);
}

static void computeLambdaVector(
    apf::Mesh* mesh,
    apf::MeshEntity* e,
    apf::Field* f,
    apf::Field* fp1,
    apf::Field* THETA_Field,
    mth::Vector<double>& lambda)
{
  apf::FieldShape* fp1s = fp1->getShape();
  int order = fp1s->getOrder();
  int etype = mesh->getType(e);
  PCU_ALWAYS_ASSERT(etype == apf::Mesh::TET);
  int nedofs = apf::countElementNodes(fp1s, etype);
  lambda.resize(nedofs);


  // get the downward faces of the element
  apf::Downward faces;
  int nf = mesh->getDownward(e, 2, faces);

  lambda.zero();
  // assemble lambda vector LOOP OVER DOWNWARD FACES
  for (int ii = 0; ii < nf; ii++) {
    apf::MeshEntity* face = faces[ii];

    // 1. get upward tets of the current face
    apf::Up up;
    mesh->getUp(face, up);
    if (crv::isBoundaryEntity(mesh, face))
      PCU_ALWAYS_ASSERT( up.n == 1);
    else
      PCU_ALWAYS_ASSERT( up.n == 2);

    apf::MeshEntity* firstTet = up.e[0];
    apf::MeshEntity* secondTet;
    if (up.n == 2)
      secondTet  = up.e[1];

    // 2. get downward edges of the face
    apf::Downward edges;
    int nedges = mesh->getDownward(face, 1, edges);

    // 3. get theta coeffs on the face
    double components[3];
    apf::getComponents(THETA_Field, face, 0, components);
    mth::Vector<double> theta_coeffs(3);
    theta_coeffs(0) = components[0];
    theta_coeffs(1) = components[1];
    theta_coeffs(2) = components[2];

    int ftype = mesh->getType(face);
    PCU_ALWAYS_ASSERT(ftype == apf::Mesh::TRIANGLE);
    int nfdofs = apf::countElementNodes(f->getShape(), ftype);
    apf::NewArray<apf::Vector3> vectorshape(nfdofs);

    apf::MeshElement* fme = apf::createMeshElement(mesh, face);
    apf::Element* fel = apf::createElement(f, fme);
    int np = apf::countIntPoints(fme, 2*order);

    // 4. Compute integral on the face
    apf::Vector3 p, tet1xi, tet2xi, curl1, curl2, curl,
      fnormal1, fnormal2, tk, vshape;
    for (int n = 0; n < np; n++) {

      apf::getIntPoint(fme, 2*order, n, p);
      double weight = apf::getIntWeight(fme, 2*order, n);
      apf::Matrix3x3 fJ;
      apf::getJacobian(fme, p, fJ);
      double jdet = apf::getJacobianDeterminant(
          fJ, apf::getDimension(mesh, face));

      // evaluate theta vector using theta coeffs
      apf::Vector3 theta_vector;
      theta_vector.zero();
      apf::NewArray<apf::Vector3> triVectorShapes (nfdofs);
      apf::getVectorShapeValues(fel, p, triVectorShapes);
      for (int i = 0; i < nedges; i++) {
        apf::Vector3 v = triVectorShapes[i];
        v = v * theta_coeffs[i];
        theta_vector += v;
      }

      // compute face outward normals wrt tets
      if (e == firstTet)
        fnormal1 = computeFaceOutwardNormal(mesh, firstTet, face, p);
      else
        fnormal1 = computeFaceOutwardNormal(mesh, secondTet, face, p);
      if (up.n == 2) {
        if (e == firstTet)
          fnormal2 = computeFaceOutwardNormal(mesh, secondTet, face, p);
        else {
          fnormal2 = computeFaceOutwardNormal(mesh, firstTet, face, p);
          theta_vector = theta_vector * -1.;
        }
      }

      curl.zero();
      // compute curl1
      tet1xi = apf::boundaryToElementXi(mesh, face, firstTet, p);
      apf::MeshElement* me1 = apf::createMeshElement(mesh, firstTet);
      apf::Element* el1 = apf::createElement(f, me1);
      apf::getCurl(el1, tet1xi, curl1);
      apf::Vector3 temp1 = apf::cross(fnormal1, curl1);
      curl += temp1;
      apf::destroyElement(el1);
      apf::destroyMeshElement(me1);

      // compute curl2
      if (up.n == 2) {
        tet2xi = apf::boundaryToElementXi(mesh, face, secondTet, p);
        apf::MeshElement* me2 = apf::createMeshElement(mesh, secondTet);
        apf::Element* el2 = apf::createElement(f, me2);
        apf::getCurl(el2, tet2xi, curl2);
        apf::Vector3 temp2 = apf::cross(fnormal2, curl2);
        curl += (temp2 * -1.);
        curl = curl * 1./2.;
        apf::destroyElement(el2);
        apf::destroyMeshElement(me2);
      }

      // compute tk (inter-element averaged flux)
      tk = curl;

      // compute p+1 order 3D vector shapes
      apf::NewArray<apf::Vector3> tetVectorShapes (nedofs);
      apf::MeshElement* me = apf::createMeshElement(mesh, e);
      apf::Element* el = apf::createElement(fp1, me);
      apf::Vector3 tetxi = apf::boundaryToElementXi(mesh, face, e, p);
      apf::getVectorShapeValues(el, tetxi, tetVectorShapes);

      apf::destroyElement(el);
      apf::destroyMeshElement(me);

      // compute integral
      apf::Vector3 theta_plus_tk = theta_vector + tk;
      double w = weight * jdet;
      theta_plus_tk = theta_plus_tk * w;

      // matrix vector multiplication
      for (int i = 0; i < nedofs; i++)
        lambda(i) += theta_plus_tk * tetVectorShapes[i];

    } // end integral loop
    apf::destroyMeshElement(fme);
  } // end face loop
}

static void getEssentialElementNDDofs(apf::Mesh* mesh, apf::MeshEntity* e,
    apf::Field* f, apf::NewArray<int>& ess_dofs, apf::NewArray<int>& uness_dofs)
{
  apf::FieldShape* fs = f->getShape();
  int etype = mesh->getType(e);
  PCU_ALWAYS_ASSERT(etype == apf::Mesh::TET);
  int nedofs = apf::countElementNodes(fs, etype);

  apf::NewArray<int> marker_list (nedofs);
  for (int i = 0; i < nedofs; i++)
    marker_list[i] = 0;

  // populate marker list by iterating over downward edges and faces
  apf::Downward edges;
  int nedges = mesh->getDownward(e, 1, edges);
  PCU_ALWAYS_ASSERT(nedges == 6);
  for (int i = 0; i < nedges; i++) {
    int nodesOnEdge = fs->countNodesOn(mesh->getType(edges[i]));
    if ( crv::isBoundaryEntity(mesh, edges[i]) ) {
      for (int n = 0; n < nodesOnEdge; n++) {
        marker_list[(i*nodesOnEdge)+n] = -1;
      }
    }
  }
  int nodesOnEdges = fs->countNodesOn(apf::Mesh::EDGE) * nedges;

  apf::Downward faces;
  int nfaces = mesh->getDownward(e, 2, faces);
  PCU_ALWAYS_ASSERT(nfaces == 4);
  for (int i = 0; i < nfaces; i++) {
    int nodesOnFace = fs->countNodesOn(mesh->getType(faces[i]));
    if ( crv::isBoundaryEntity(mesh, faces[i]) ) {
      for (int n = 0; n < nodesOnFace; n++) {
        marker_list[(nodesOnEdges + i*nodesOnFace)+n] = -1;
      }
    }
  }

  int num_marked = 0;
  for (int i = 0; i < nedofs; i++) {
    if (marker_list[i])
      num_marked++;
  }

  // use marker list to get ess_dofs list
  ess_dofs.allocated() ? ess_dofs.resize(num_marked)
    : ess_dofs.allocate(num_marked);
  uness_dofs.allocated() ? uness_dofs.resize(nedofs-num_marked)
    : uness_dofs.allocate(nedofs-num_marked);
  int ess_dof_counter = 0;
  int uness_dof_counter = 0;
  for (int i = 0; i < nedofs; i++) {
    if(marker_list[i])
      ess_dofs[ess_dof_counter++] = i;
    else
      uness_dofs[uness_dof_counter++] = i;
  }
}

/**
 * Inputs: Matrix A, Vector X, Vector B, essential dofs, other dofs.
 * Output: reduced matrix A, reduced rhs B.
 */
static void eliminateDBCs(
    mth::Matrix<double> const &A,
    mth::Vector<double> const &X,
    mth::Vector<double> const &B,
    apf::NewArray<int>  const &ess_dofs,
    apf::NewArray<int>  const &uness_dofs,
    mth::Matrix<double> &Anew,
    mth::Vector<double> &Bnew)
{
  int num_ess_dofs = ess_dofs.size();
  int num_uness_dofs = uness_dofs.size();

  // 1. Remove rows of A corresponding to
	// ess_dofs by copying into Anew
  Anew.resize(num_uness_dofs, num_uness_dofs);
	for(int rr = 0; rr < num_uness_dofs; rr++) {
    int i = uness_dofs[rr];
    for(int cc = 0; cc < num_uness_dofs; cc++) {
      int j = uness_dofs[cc];
      Anew(rr,cc) = A(i,j);
    }
  }

	// 2. Assemble new B
	Bnew.resize(num_uness_dofs);
  for(int i = 0; i < num_uness_dofs; i++) {
    Bnew(i) = B(uness_dofs[i]);
  }

  if (num_ess_dofs > 0) {
	  // 3. Subtract from Bnew: (Bnew -= Ae*Xe)
    mth::Matrix<double> Ae(num_uness_dofs, num_ess_dofs);
    for(int rr = 0; rr < num_uness_dofs; rr++) {
      int i = uness_dofs[rr];
      for(int cc = 0; cc < num_ess_dofs; cc++) {
        int j = ess_dofs[cc];
        Ae(rr,cc) = A(i,j);
      }
    }

    mth::Vector<double> Xe(num_ess_dofs);
    for(int i = 0; i < num_ess_dofs; i++) {
      Xe(i) = X(ess_dofs[i]);
    }

    mth::Vector<double> temp;
    mth::multiply(Ae, Xe, temp);
    Bnew -= temp;
  }
}

static double computeL2Error(apf::Mesh* mesh, apf::MeshEntity* e,
	apf::Field* f, mth::Vector<double> const error_dofs)
{
	double error = 0.0;

  apf::FieldShape* fs = f->getShape();
  int type = mesh->getType(e);
	PCU_ALWAYS_ASSERT(type == apf::Mesh::TET);
  int nd = apf::countElementNodes(fs, type);
  int dim = apf::getDimension(mesh, e);
  double w;

  apf::MeshElement* me = apf::createMeshElement(mesh, e);
  apf::Element* el = apf::createElement(f, me);
  int order = 2 * fs->getOrder();
  int np = apf::countIntPoints(me, order);

  apf::NewArray<apf::Vector3> vectorshape(nd);

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
    for (int j = 0; j < nd; j++)
      for (int k = 0; k < dim; k++)
        vectorShape(j,k) = vectorshape[j][k];

    mth::Matrix<double> vectorShapeT (dim, nd);
    mth::transpose(vectorShape, vectorShapeT);

		mth::Vector<double> err_func;
		mth::multiply(vectorShapeT, error_dofs, err_func);

		error += w * (err_func * err_func);
	}
	if (error < 0.0)
		error = -error;

	return sqrt(error);
}

apf::Field* computeErrorField(apf::Field* ef, apf::Field* THETA_Field)
{

	// 1. Create per-element SCALAR error field
  apf::Field* error_field = apf::createIPField(
		apf::getMesh(ef), "residual_error_field", apf::SCALAR, 1);

  // 2. Create p+1 order tet ND field
  int order = ef->getShape()->getOrder();
  int orderp1 = order+1;
  apf::Field* efp1 = apf::createField(apf::getMesh(ef),
	 "orderp1_nedelec_field", apf::SCALAR, apf::getNedelec(orderp1));
  apf::zeroField(efp1);

  // 2. iterate over all elements of the mesh
  apf::MeshEntity* el;
  apf::MeshIterator* it = apf::getMesh(ef)->begin(3);
  while ((el = apf::getMesh(ef)->iterate(it))) {

    // 2(a). Assemble LHS element matrix
    mth::Matrix<double> A;
    assembleElementMatrix( apf::getMesh(ef), el, efp1, A);

    // 2(b). Compute Bilinear Form Vector
    mth::Vector<double> blf;
    computeResidualBLF(apf::getMesh(efp1), el, ef, efp1, blf);

    // 2(c). Compute Linear Form Vector
    mth::Vector<double> lf;
    assembleDomainLFElementVector(apf::getMesh(efp1), el, efp1, lf);

    // 2(d). Compute Lambda Vector
    mth::Vector<double> lambda;
    computeLambdaVector(
        apf::getMesh(ef), el, ef, efp1, THETA_Field, lambda);

    // 2(e). Assemble RHS element vector = blf - lf - lambda
    mth::Vector<double> B(blf.size());
    B.zero();
    B += blf; B -= lf; B -= lambda;

    // 2(f). Get List of Essential Dofs
    apf::NewArray<int> ess_dofs, uness_dofs;
    getEssentialElementNDDofs(
        apf::getMesh(efp1), el, efp1, ess_dofs, uness_dofs);

    // 2(g). eliminate Dirichlet (Essential) Boundary Conditions
		mth::Vector<double> X, Bnew;
    mth::Matrix<double> Anew;
		X.resize(B.size());
		X.zero(); // initialize X with exact DBC (e = 0.0)
		eliminateDBCs(A, X, B, ess_dofs, uness_dofs, Anew, Bnew);

		// 2(h). Solve the reduced system
	  mth::Matrix<double> Q, R;
		mth::decomposeQR(Anew, Q, R);
		mth::Vector<double> Xnew;
		mth::solveFromQR(Q, R, Bnew, Xnew);

		// 2(i). Recover the solution
		mth::Vector<double> error_dofs(B.size());
    for(unsigned int i = 0; i < ess_dofs.size(); i++) {
    	int index = ess_dofs[i];
      error_dofs(index) = X(index);
    }
    for(unsigned int i = 0; i < uness_dofs.size(); i++) {
      int index = uness_dofs[i];
      error_dofs(index) = Xnew(i);
    }

		// 2(j). Compute L2 Norm Error
		double l2_error = computeL2Error(apf::getMesh(ef), el, efp1, error_dofs);
		apf::setScalar(error_field, el, 0, l2_error);
  }
  apf::getMesh(ef)->end(it);
	apf::destroyField(efp1);

	return error_field;
}

apf::Field* estimateError(apf::Field* f)
{
  apf::Field* g = em::equilibrateResiduals(f);
  apf::Field* THETA = em::computeFluxCorrection(f, g);
  apf::destroyField(g);

  apf::Field* error_field = em::computeErrorField(f, THETA);
  apf::destroyField(THETA);

  return error_field;
}


/*class SizeFieldOp : public apf::CavityOp
{
public:
  SizeFieldOpOp(apf::Field* error_field):
    apf::CavityOp(error_field->getMesh())
  {}
  virtual Outcome setEntity(apf::MeshEntity* e)
  {
    entity = e;
    if (apf::hasEntity(error_field, entity))
      return SKIP;
    if ( !requestLocality(&entity,1))
      return REQUEST;
    return OK;
  }
  virtual void apply()
  {

  }
  apf::MeshEntity* e;
};*/

}

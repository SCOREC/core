/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include <iostream>
#include <cstdlib>
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

#include "em.h"

namespace em {
//TODO destroy all fields created inside functions to prevent memory leaks.

static void computeResidualBLF(apf::Mesh* mesh, apf::MeshEntity* e,
  apf::Field* f, mth::Vector<double>& blf)
{
  apf::FieldShape* fs = f->getShape();
  int type = mesh->getType(e);
  PCU_ALWAYS_ASSERT(type == apf::Mesh::TET);
  int nd = apf::countElementNodes(fs, type);
  int dim = apf::getDimension(mesh, e);
  double w;

  apf::NewArray<apf::Vector3> curlshape(nd);
  apf::NewArray<apf::Vector3> vectorshape(nd);
  mth::Matrix<double> phys_curlshape(nd, dim);
  mth::Matrix<double> Vectorshape(nd, dim);
  mth::Vector<double> curlcurl_vec (nd);
  mth::Vector<double> mass_vec (nd);
  blf.resize(nd);

  apf::MeshElement* me = apf::createMeshElement(mesh, e);
  apf::Element* el = apf::createElement(f, me);
  int order = 2 * fs->getOrder() - 2;
  int np = apf::countIntPoints(me, order); // int points required

  // 1. Compute Curl Curl Integration
  curlcurl_vec.zero();
  apf::Vector3 p, curl;
  for (int i = 0; i < np; i++) {
    apf::getIntPoint(me, order, i, p);
    double weight = apf::getIntWeight(me, order, i);
    apf::Matrix3x3 J;
    apf::getJacobian(me, p, J);
    double jdet = apf::getJacobianDeterminant(J, dim);
    w = weight; // TODO check why do not need division by jdet

    // get curl vector
    apf::getCurl(el, p, curl);
    
    // get curlshape values
    el->getShape()->getLocalVectorCurls(mesh, e, p, curlshape);
    phys_curlshape.zero();
    for (int i = 0; i < nd; i++)
      for (int j = 0; j < dim; j++)
        for (int k = 0; k < dim; k++)
          phys_curlshape(i,j) += curlshape[i][k] * J[k][j];

    // multiply
    mth::Vector<double> temp (nd);
    temp.zero();
    mth::Vector<double> c (dim);
    c(0) = curl[0]; c(1) = curl[1]; c(2) = curl[2];
    mth::multiply(phys_curlshape, c, temp);
    temp *= w;
    
    curlcurl_vec += temp;

  }
    
  // 2. Compute Vector Mass Integration
  mass_vec.zero();
  apf::Vector3 vvalue;
  for (int i = 0; i < np; i++) {
    apf::getIntPoint(me, order, i, p);
    double weight = apf::getIntWeight(me, order, i);
    apf::Matrix3x3 J;
    apf::getJacobian(me, p, J);
    double jdet = apf::getJacobianDeterminant(J, dim);

    w = weight * jdet;

    apf::getVector(el, p, vvalue);

    apf::getVectorShapeValues(el, p, vectorshape);
    Vectorshape.zero();
    for (int i = 0; i < nd; i++)
      for (int j = 0; j < dim; j++)
        Vectorshape(i,j) = vectorshape[i][j];


    mth::Vector<double> temp (nd);
    temp.zero();
    mth::Vector<double> v (dim);
    v(0) = vvalue[0]; v(1) = vvalue[1]; v(2) = vvalue[2];
    mth::multiply(Vectorshape, v, temp);
    temp *= w;

    mass_vec += temp;
  }

  // 3. get result
  blf.zero();
  blf += curlcurl_vec;
  blf += mass_vec;

  // TODO Take care of Negative Dofs 
}

// 
// @ apf::Field* f --> input p-order ND field
// @ apf::Field* fp1 --> p+1 order ND field created for local BVPs
// @ apf::Field* THETA --> constant face field containing theta coeffs
static void computeLambdaVector(apf::Mesh* mesh, apf::MeshEntity* e,
  apf::Field* f, apf::Field* fp1, apf::Field* THETA, mth::Vector<double>& lambda)
{
  apf::FieldShape* fp1s = fp1->getShape();
  int etype = mesh->getType(e);
  PCU_ALWAYS_ASSERT(etype == apf::Mesh::TET);
  int nedofs = apf::countElementNodes(fp1s, etype);
  int edim = apf::getDimension(mesh, e);
  lambda.resize(nedofs);

  // create a 2D nedelec field 
  int order = fp1s->getOrder();
  apf::Field* faceNDField = apf::createField(
      mesh, "face_nedelec_field", apf::SCALAR, apf::getNedelec(order));
  apf::zeroField(faceNDField);
  apf::FieldShape* faceNDFieldShape = faceNDField->getShape();

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
    apf::getComponents(THETA, face, 0, components);
    mth::Vector<double> theta_coeffs(components); // TODO Warning make this work

    int ftype = mesh->getType(face);
    PCU_ALWAYS_ASSERT(ftype == apf::Mesh::TRIANGLE);
    int nfdofs = apf::countElementNodes(faceNDFieldShape, ftype);
    int fdim = apf::getDimension(mesh, face);
    apf::NewArray<apf::Vector3> vectorshape(nfdofs);

    apf::MeshElement* fme = apf::createMeshElement(mesh, face);
    apf::Element* fel = apf::createElement(faceNDField, fme);
    int np = apf::countIntPoints(fme, 2*order); // int points required

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
      apf::NewArray<apf::Vector3> vector2Dshapes (nfdods);
      apf::getVectorShapeValues(fel, p, vector2Dshapes);
      int which, rotate; bool flip; // negative ND dofs
      for (int i = 0; i < theta_coeffs.size(); i++) {
        apf::getAlignment(mesh, face, edges[i], which, flip, rotate);
        apf::Vector3 v = vector2Dshapes[i];
        if (flip) { v = v * -1.; }

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
        else
          fnormal2 = computeFaceOutwardNormal(mesh, firstTet, face, p);
        std::cout << "normal1 " << fnormal1 << std::endl; // REMOVE
        std::cout << "normal2 " << fnormal2 << std::endl; // REMOVE
      }

      curl.zero();
      // compute curl1
      tet1xi = apf::boundaryToElementXi(mesh, face, firstTet, p);
      apf::MeshElement* me1 = apf::createMeshElement(mesh, firstTet);
      apf::Element* el1 = apf::createElement(f, me1);
      apf::getCurl(el1, tet1xi, curl1);
      apf::Vector3 temp1 = apf::cross(fnormal1, curl1); // TODO clean
      curl += temp1; //
      apf::destroyElement(el1);
      apf::destroyMeshElement(me1);

      // compute curl2
      if (up.n == 2) {
        tet2xi = apf::boundaryToElementXi(mesh, face, secondTet, p);
        apf::MeshElement* me2 = apf::createMeshElement(mesh, secondTet);
        apf::Element* el2 = apf::createElement(f, me2);
        apf::getCurl(el2, tet2xi, curl2);
        apf::Vector3 temp2 = apf::cross(fnormal2, curl2); // TODO clean
        curl += (temp2 * -1.); //
        apf::destroyElement(el2);
        apf::destroyMeshElement(me2);
      }

      // compute tk (inter-element averaged flux)
      tk = curl * 1./2.;
      std::cout << "tk " << tk << std::endl;

      // compute p+1 order 3D vector shapes
      apf::NewArray<apf::Vector3> vector3dshapes (nedofs);
      apf::MeshElement* me = apf::createMeshElement(mesh, e);
      apf::Element* el = apf::createElement(fp1, me);
      apf::Vector3 tetxi = apf::boundaryToElementXi(mesh, face, e, p);
      apf::getVectorShapeValues(el, tetxi, vector3dshapes);
      // TODO take care of negative ND dofs
      apf::destroyMeshElement(me); 
      // compute integral
      apf::Vector3 theta_plus_tk = theta_vec + tk;
      double w = weight * jdet;
      theta_plus_tk *= w;

      // matrix vector multiplication
      for (int i = 0; i < nedofs; i++)
        lambda(i) += theta_plus_tk * vector3dshapes[i];

    } // end integral loop
    apf::destroyMeshElement(fme);
  } // end face loop
}

// @ apf::Field* ef --> input p order ND electric field
// @ apf::Field* ef --> input theta corrected flux field
// *** p+1 order field is created inside this function
apf::Field* emEstimateError(apf::Field* ef, apf::Field* correctedFlux)
{
  // 1. Create one order higher ND field
  int order = ef->getShape()->getOrder();
  int orderp1 = order+1;
  apf::Field* efp1 = apf::createField(
      apf::getMesh(ef), "higher_order_nedelec_field", apf::SCALAR, apf::getNedelec(orderp1));
  apf::zeroField(efp1);

  // 2. iterate over all elements of the mesh
  apf::MeshEntity* el;
  apf::MeshIterator* it = apf::getMesh(ef)->begin(3);
  while ((el = apf::getMesh(ef)->iterate(it))) {
    // 2(a). Assemble LHS element matrix
    mth::Matrix<double> lhs;
    assembleElementMatrix( apf::getMesh(ef), el, efp1, lhs);
    // TODO Take care of negative dofs in lhs

    // 2(b). Compute Bilinear Form Vector 
    mth::Vector<double> blf;
    computeResidualBLF(apf::getMesh(efp1), e, efp1, blf);
    // TODO Take care of negative dofs in blf

    // 2(c). Compute Linear Form Vector
    mth::Vector<double> lf;
    assembleDomainLFElementVector(apf::getMesh(efp1), e, efp1, lf);
    // TODO Take care of negative dofs in lf
    
    // 2(d). Compute Lamda Vector
    mth::Vector<double> lambda;
    computeLambdaVector(apf::getMesh(efp1), e, ef, efp1, correctedFlux, lambda);
    // TODO Take care of negative dofs in lambda
   
     


  }
  apf::getMesh(ef)->end(it);



}

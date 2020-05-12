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
    
    // 2(c). Compute Lamda Vector
    mth::Vector<double> lambda;
    assembleDomainLFElementVector(apf::getMesh(efp1), e, efp1, lf);
    // TODO Take care of negative dofs in lambda


  }
  apf::getMesh(ef)->end(it);



}

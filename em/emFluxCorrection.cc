/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include <iostream>
#include <cstdlib>

#include "em.h"

namespace em {

struct QRDecomp {
  mth::Matrix<double> Q;
  mth::Matrix<double> R;
};

static void assembleFaceMassMatrix(apf::Mesh* mesh,apf::MeshEntity* e,
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

  apf::Downward edges;
  int ned =  mesh->getDownward(e, 1, edges);
  PCU_ALWAYS_ASSERT(ned == 3);
  int which, rotate; bool flip;

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
    for (int j = 0; j < nd; j++)
      for (int k = 0; k < dim; k++)
        vectorShape(j,k) = vectorshape[j][k];

    // negate negatve dof indices
    for (int ei = 0; ei < ned; ei++) {
      apf::getAlignment(mesh, e, edges[ei], which, flip, rotate);
      if(flip) {
        for (int j = 0; j < dim; j++)
          vectorShape(ei, j) = -1*vectorShape(ei, j);
      }
    }

    mth::Matrix<double> vectorShapeT (dim, nd);
    mth::transpose(vectorShape, vectorShapeT);
    mth::Matrix<double> M (nd,nd);
    M.zero();
    mth::multiply(vectorShape, vectorShapeT, M);
    M *= w;
    elmat += M;
  }
  apf::destroyElement(el);
  apf::destroyMeshElement(me);
}


// This function takes the NedelecField and computes
// correction vectors on faces and stores them in a field
apf::Field* computeFluxCorrection(apf::Field* ef, apf::Field* g)
{
  int order = ef->getShape()->getOrder();
  apf::Field* faceNedelecField = apf::createField(
      apf::getMesh(g), "face_nedelec_field", apf::SCALAR, apf::getNedelec(order));
  apf::zeroField(faceNedelecField);


  apf::Field* correctedFluxField =  createPackedField(
      apf::getMesh(g), "corrected_flux_field", 3, apf::getConstant(2));

  // iterate over all faces of the mesh
  apf::MeshEntity* face;
  apf::MeshIterator* it = apf::getMesh(g)->begin(2);
  while ((face = apf::getMesh(g)->iterate(it))) {
    // assemble RHS vector
    double components[3];
    apf::getComponents(g, face, 0, components);
    mth::Vector<double> rhs(3);
    rhs(0) = components[0]; rhs(1) = components[1]; rhs(2) = components[2];

    // assemlbe face mass matrix
    mth::Matrix<double> M;
    assembleFaceMassMatrix(
        apf::getMesh(g), face, faceNedelecField, M);

    // solve the system
    QRDecomp qr;
    mth::decomposeQR(M, qr.Q, qr.R);
    std::cout << "M" << std::endl;
    std::cout << M << std::endl;
    std::cout << "Q" << std::endl;
    std::cout << qr.Q << std::endl;
    std::cout << "R" << std::endl;
    std::cout << qr.R << std::endl;
    std::cout << "RHS" << std::endl;
    std::cout << rhs << std::endl;
    
    mth::Vector<double> theta;
    mth::solveFromQR(qr.Q, qr.R, rhs, theta);
    std::cout << "theta" << std::endl;
    std::cout << theta << std::endl;

    // set solution vector on face field
    //theta.toArray(components);
    components[0] = theta(0);
    components[1] = theta(1);
    components[2] = theta(2);
    apf::setComponents(
        correctedFluxField, face, 0, components);
  }
  apf::getMesh(g)->end(it);

  return correctedFluxField;
}

}



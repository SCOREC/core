/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include <iostream>
#include <cstdlib>

#include "crv.h"
#include "crvShape.h"
#include "apfElement.h"

#include "em.h"
using namespace std;
namespace em {

struct QRDecomp {
  mth::Matrix<double> Q;
  mth::Matrix<double> R;
};

static void assembleFaceMassMatrix(apf::Mesh* mesh, apf::MeshEntity* e,
    apf::Field* f, mth::Matrix<double>& elmat)
{
  apf::FieldShape* fs = f->getShape();
  int type = mesh->getType(e);
  int nd = apf::countElementNodes(fs, type);
  int dim = apf::getDimension(mesh, e);
  int sdim = mesh->getDimension();
  PCU_ALWAYS_ASSERT(type == apf::Mesh::TRIANGLE && dim == 2);
  double w;

  apf::NewArray<apf::Vector3> vectorshapes(nd);
  elmat.resize(nd,nd);

  apf::MeshElement* me = apf::createMeshElement(mesh, e);
  apf::Element* el = apf::createElement(f, me);
  int int_order = 2 * fs->getOrder();
  int np = apf::countIntPoints(me, int_order); // int points required

  apf::Downward edges;
  int ned =  mesh->getDownward(e, 1, edges);
  int which, rotate; bool flip;

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

    // negate negative dof indices
    // TODO maybe can do this outside the function and therfore
    // merge this function with assembleVectorMassElementMatrix
    for (int ei = 0; ei < ned; ei++) {
      apf::getAlignment(mesh, e, edges[ei], which, flip, rotate);
      if(flip) {
        for (int j = 0; j < sdim; j++)
          vectorShapes(ei, j) = -1*vectorShapes(ei, j);
      }
    }

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


// This function takes the NedelecField and computes
// correction vectors on faces and stores them in a field
apf::Field* computeFluxCorrection(apf::Field* ef, apf::Field* g)
{
  int order = ef->getShape()->getOrder();
  cout << "tri_field_order " << order << endl; // REMOVE
  apf::Field* triNedelecField = apf::createField(
      apf::getMesh(ef), "tri_nedelec_field", apf::SCALAR, apf::getNedelec(order));
  apf::zeroField(triNedelecField);

  apf::Field* correctedFluxField =  createPackedField(
      apf::getMesh(ef), "corrected_flux_field", 3, apf::getConstant(2));

  // iterate over all faces of the mesh
  apf::MeshEntity* face;
  apf::MeshIterator* it = apf::getMesh(ef)->begin(2);
  while ((face = apf::getMesh(ef)->iterate(it))) {
    // 1. assemble RHS vector
    double components[3];
    apf::getComponents(g, face, 0, components);
    mth::Vector<double> rhs(3);
    rhs(0) = components[0]; rhs(1) = components[1]; rhs(2) = components[2];

    bool debug = true;
    if (debug) {
      cout << "RHS VECTOR" << endl;
      cout << components[0] << " " << components[1] << " " << components[2] << endl;
      cout << rhs << endl;
      cout << "==============" << endl;
    }

    // 2. assemble face mass matrix
    mth::Matrix<double> M;
    assembleFaceMassMatrix(
        apf::getMesh(ef), face, triNedelecField, M);

    // 3. solve the system
    QRDecomp qr;
    mth::decomposeQR(M, qr.Q, qr.R);
    std::cout << "M" << std::endl; // REMOVE DEBUG
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
    //theta.toArray(components); TODO clean
    components[0] = theta(0);
    components[1] = theta(1);
    components[2] = theta(2);
    apf::setComponents(
        correctedFluxField, face, 0, components);
  }
  apf::getMesh(ef)->end(it);

  apf::destroyField(triNedelecField);

  return correctedFluxField;
}

}



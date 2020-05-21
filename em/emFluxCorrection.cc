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


// This function takes the NedelecField and computes
// correction vectors on faces and stores them in a field
apf::Field* computeFluxCorrection(apf::Field* ef, apf::Field* g)
{
  apf::Mesh* mesh = apf::getMesh(ef);
  int order = ef->getShape()->getOrder();

  apf::Field* triNedelecField = apf::createField(
      mesh, "tri_nedelec_field", apf::SCALAR, apf::getNedelec(order));
  apf::zeroField(triNedelecField);

  apf::Field* correctedFluxField =  createPackedField(
      mesh, "corrected_flux_field", 3, apf::getConstant(2));

  apf::MeshEntity* face;
  apf::MeshIterator* it = mesh->begin(2);
  while ((face = mesh->iterate(it))) {

    // 1. assemble RHS vector
    double components[3];
    apf::getComponents(g, face, 0, components);
    cout << "initial gs without negation" << endl; // REMOVE
    cout << components[0] << " " << components[1] << " " << components[2] << endl; // REMOVE

    apf::Downward edges;
    int ne = mesh->getDownward(face, 1, edges);
    for (int i = 0; i < ne; i++) {
      apf::setScalar(triNedelecField, edges[i], 0, components[i]);
    }

    apf::MeshElement* me = apf::createMeshElement(mesh, face);
    apf::Element* el = apf::createElement(triNedelecField, me);
    apf::NewArray<double> facegs;
    el->getElementDofs(facegs);

    mth::Vector<double> rhs(facegs.size()); // TODO clean
    for (size_t i = 0; i < facegs.size(); i++) {
      rhs(i) = facegs[i];
    }
    cout << "final gs with negation" << endl; // REMOVE
    cout << rhs << endl; // REMOVE

    apf::destroyElement(el);
    apf::destroyMeshElement(me);

    // 2. assemble LHS face mass matrix
    mth::Matrix<double> M;
    assembleFaceMassMatrix(
        mesh, face, triNedelecField, M);

    // 3. solve the system
    QRDecomp qr;
    mth::decomposeQR(M, qr.Q, qr.R);
    mth::Vector<double> theta;
    mth::solveFromQR(qr.Q, qr.R, rhs, theta);
    std::cout << "LHS Face Mass Matrix" << std::endl; // REMOVE DEBUG
    std::cout << M << std::endl;
    std::cout << "RHS Vector" << std::endl;
    std::cout << rhs << std::endl;
    std::cout << "theta coeffs" << std::endl;
    std::cout << theta << std::endl;

    // set solution vector on face field
    //theta.toArray(components); TODO clean
    components[0] = theta(0);
    components[1] = theta(1);
    components[2] = theta(2);
    apf::setComponents(
        correctedFluxField, face, 0, components);
  }
  mesh->end(it);

  apf::destroyField(triNedelecField);

  return correctedFluxField;
}

}



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

/*
 * This function takes the NedelecField and g parameters field
 * and returns the THETA field. It computes the 3 scalar parameters
 * per face for the flux correction vectors.
 */
apf::Field* computeFluxCorrection(apf::Field* ef, apf::Field* g)
{
  apf::Mesh* mesh = apf::getMesh(ef);
  apf::Field* THETA_Field =  createPackedField(
      mesh, "theta_field", 3, apf::getConstant(2));

  apf::MeshEntity* face;
  apf::MeshIterator* it = mesh->begin(2);
  while ((face = mesh->iterate(it))) {

    // 1. assemble RHS vector
    double components[3];
    apf::getComponents(g, face, 0, components);

    mth::Vector<double> rhs(3);
    apf::Downward edges;
    int ne = mesh->getDownward(face, 1, edges);
    for (int i = 0; i < ne; i++) {
      int which, rotate; bool flip;
      apf::getAlignment(mesh, face, edges[i], which, flip, rotate);
      if (flip)
        rhs(i) = -1. * components[i];
      else
        rhs(i) = components[i];
    }

    // 2. assemble LHS face mass matrix
    mth::Matrix<double> M;
    assembleVectorMassElementMatrix(mesh, face, ef, M);

    // 3. solve the system
    QRDecomp qr;
    mth::decomposeQR(M, qr.Q, qr.R);
    mth::Vector<double> theta;
    mth::solveFromQR(qr.Q, qr.R, rhs, theta);

    // set solution vector on face field
    components[0] = theta(0);
    components[1] = theta(1);
    components[2] = theta(2);
    apf::setComponents(THETA_Field, face, 0, components);
  }
  mesh->end(it);

  return THETA_Field;
}

}



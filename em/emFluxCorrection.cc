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

struct QRDecomp {
  mth::Matrix<double> Q;
  mth::Matrix<double> R;
};

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
    assembleVectorMassElementMatrix(
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



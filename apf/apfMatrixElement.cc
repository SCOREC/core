/*
 * Copyright 2018 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfMatrixField.h"
#include "apfMatrixElement.h"
#include "apfVectorElement.h"

namespace apf {

MatrixElement::MatrixElement(MatrixField* f, MeshElement* e):
    ElementOf<Matrix3x3>(f,e)
{
}

MatrixElement::~MatrixElement()
{
}

// laid out in array as F_i*3+j+9*d
void MatrixElement::grad(Vector3 const& xi, Vector<27>& g)
{
  Matrix3x3* nodeValues = getNodeValues();
  NewArray<Vector3> globalGradients;
  getGlobalGradients(xi, globalGradients);
  // for the first time through g, set the values of g
  for(int i=0; i<3; ++i) {
    for(int j=0; j<3; ++j) {
      for(int d=0; d<3; ++d) {
          g[i*3+j+d*9]= nodeValues[0][i][j]*globalGradients[0][d];
      }
    }
  }
  for(int nd=1; nd<nen; ++nd) {
    for(int i=0; i<3; ++i) {
      for(int j=0; j<3; ++j) {
        for(int d=0; d<3; ++d) {
            g[i*3+j+d*9]+= nodeValues[nd][i][j]*globalGradients[nd][d];
        }
      }
    }
  }
}
}//namespace apf

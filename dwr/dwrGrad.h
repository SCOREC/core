/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef DWR_GRAD_H
#define DWR_GRAD_H

namespace dwr {

template<class Matrix3x3T, class Vector3T>
void computeVectorGrad(
    int numDims,
    int numNodes,
    apf::DynamicArray<Vector3T> const& u,
    apf::NewArray<apf::Vector3> const& gradBF,
    Matrix3x3T& gradU)
{
  zeroMatrix3x3(gradU);
  for (int i=0; i < numNodes; ++i)
  for (int c=0; c < numDims; ++c)
  for (int d=0; d < numDims; ++d)
    gradU[c][d] += u[i][c] * gradBF[i][d];
}

}

#endif

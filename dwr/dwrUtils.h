/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef DWR_UTILS_H
#define DWR_UTILS_H

namespace apf {
class Mesh;
class MeshEntity;
}

namespace dwr {

template<class VectorT>
void zeroVector(int n, VectorT& v)
{
  for (int i=0; i < n; ++i)
    v[i] = 0.0;
}

template<class Matrix3x3T>
void zeroMatrix3x3(Matrix3x3T& m)
{
  for (int i=0; i < 3; ++i)
  for (int j=0; j < 3; ++j)
    m[i][j] = 0.0;
}

double getMeshSize(apf::Mesh* m, apf::MeshEntity* e);

void print(const char* format, ...) __attribute__((format(printf,1,2)));

}

#endif

/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfElement.h"
#include "apfShape.h"
#include "apfMesh.h"
#include "apfComplexType.h"
#include "apfVectorElement.h"

namespace apf {

template <class T>
ElementBase<T>::ElementBase(FieldBase* f, MeshEntity* e)
{
  init(f,e,0);
}

template <class T>
ElementBase<T>::ElementBase(FieldBase* f, VectorElement* p)
{
  init(f,p->getEntity(),p);
}

template <class T>
void ElementBase<T>::getGlobalGradients(Vector3 const& local, NewArray<Vector3>& globalGradients)
{
  Matrix3x3 J;
  parent->getJacobian(local,J);
  Matrix3x3 jinv = getJacobianInverse(J, getDimension());
  NewArray<Vector3> localGradients;
  shape->getLocalGradients(mesh, entity, local,localGradients);
  globalGradients.allocate(nen);
  for (int i=0; i < nen; ++i)
    globalGradients[i] = jinv * localGradients[i];
}

template class ElementBase<double>;
template class ElementBase<double_complex>;

Matrix3x3 getJacobianInverse(Matrix3x3 J, int dim)
{
  switch (dim) {
/* this routine computes the Moore-Penrose pseudo-inverse
   of J. We need it to handle non-square Jacobians that
   arise from edges embedded in 2D or higher and faces
   embedded in 3D. The resulting inverse Jacobian should
   be adequate for computing global shape function
   and field gradients based on local gradients */
    case 1: {
      Matrix3x3 jinvt;
      jinvt[0] = J[0] / (J[0] * J[0]);
      jinvt[1] = jinvt[2] = Vector3(0,0,0);
      return transpose(jinvt);
    }
    case 2: {
      Matrix<2,3> A;
      A[0] = J[0];
      A[1] = J[1];
      Matrix<3,2> At = transpose(A);
      Matrix<2,3> Ainvt = transpose(At * invert(A * At));
      Matrix3x3 jinvt;
      jinvt[0] = Ainvt[0];
      jinvt[1] = Ainvt[1];
      jinvt[2] = Vector3(0,0,0);
      return transpose(jinvt);
    }
    case 3:
      return invert(J);
    default:
      fail("getJacobianInverse: bad dimension");
  }
}

}//namespace apf

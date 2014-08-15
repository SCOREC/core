/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfVectorElement.h"
#include "apfVectorField.h"

namespace apf {

VectorElement::VectorElement(VectorField* f, MeshEntity* e):
  ElementOf<Vector3>(f,e)
{
}

VectorElement::VectorElement(VectorField* f, VectorElement* p):
  ElementOf<Vector3>(f,p)
{
}

double VectorElement::div(Vector3 const& xi)
{
  NewArray<Vector3> globalGradients;
  getGlobalGradients(xi,globalGradients);
  Vector3* nodeValues = getNodeValues();
  double d = globalGradients[0] * nodeValues[0];
  for (int i=1; i < nen; ++i)
    d = d + globalGradients[i] * nodeValues[i];
  return d;
}

void VectorElement::curl(Vector3 const& xi, Vector3& c)
{
  NewArray<Vector3> globalGradients;
  getGlobalGradients(xi,globalGradients);
  Vector3* nodeValues = getNodeValues();
  c = cross(globalGradients[0],nodeValues[0]);
  for (int i=1; i < nen; ++i)
    c = c + cross(globalGradients[i],nodeValues[i]);
}

void VectorElement::gradHelper(
    NewArray<Vector3>& nodalGradients,
    Matrix3x3& g)
{
  Vector3* nodeValues = getNodeValues();
  g = tensorProduct(nodalGradients[0],nodeValues[0]);
  for (int i=1; i < nen; ++i)
    g = g + tensorProduct(nodalGradients[i],nodeValues[i]);
}

void VectorElement::grad(Vector3 const& xi, Matrix3x3& g)
{
  NewArray<Vector3> globalGradients;
  getGlobalGradients(xi,globalGradients);
  gradHelper(globalGradients,g);
}

void VectorElement::getJacobian(Vector3 const& xi, Matrix3x3& J)
{
  NewArray<Vector3> localGradients;
  this->shape->getLocalGradients(xi,localGradients);
  gradHelper(localGradients,J);
}

double getJacobianDeterminant(Matrix3x3 const& J, int dimension)
{
  if (dimension == 3)
  {
    /* det(J) is also the triple product of the
       "tangent vectors" in 3D, the volume of their
       parallelpiped, which is the differential volume
       of the coordinate field */
    return getDeterminant(J);
  }
  if (dimension == 2)
  {
    /* |\frac{\partial x}{\partial s}\times
        \frac{\partial x}{\partial t}|,
        the area spanned by the tangent vectors
        at this point, surface integral. */
    return cross(J[0],J[1]).getLength();
  }
  // assuming at this point dimension=1
  /* \|\vec{x}_{,\xi}\| the length
   of the tangent vector at this point.
   line integral:
   ds = sqrt(dx^2 + dy^2 + dz^2) */
  return J[0].getLength();
}

double VectorElement::getDV(Vector3 const& xi)
{
  Matrix3x3 J;
  getJacobian(xi,J);
  return getJacobianDeterminant(J,getDimension());
}

}//namespace apf

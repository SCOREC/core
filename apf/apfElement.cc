/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfElement.h"
#include "apfShape.h"
#include "apfMesh.h"
#include "apfVectorElement.h"

namespace apf {

void Element::init(Field* f, MeshEntity* e, VectorElement* p)
{
  field = f;
  mesh = f->getMesh();
  entity = e;
  shape = f->getShape()->getEntityShape(mesh->getType(e));
  parent = p;
  nen = shape->countNodes();
  nc = f->countComponents();
  getNodeData();
}

Element::Element(Field* f, MeshEntity* e)
{
  init(f,e,0);
}

Element::Element(Field* f, VectorElement* p)
{
  init(f,p->getEntity(),p);
}

Element::~Element()
{
}

static Matrix3x3 getJacobianInverse(Matrix3x3 J, int dim)
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

void Element::getGlobalGradients(Vector3 const& local,
                                 NewArray<Vector3>& globalGradients)
{
  Matrix3x3 J;
  parent->getJacobian(local,J);
  Matrix3x3 jinv = getJacobianInverse(J, getDimension());
  NewArray<Vector3> localGradients;
  shape->getLocalGradients(local,localGradients);
  globalGradients.allocate(nen);
  for (int i=0; i < nen; ++i)
    globalGradients[i] = jinv * localGradients[i];
}

void Element::getComponents(Vector3 const& xi, double* c)
{
  NewArray<double> shapeValues;
  shape->getValues(xi,shapeValues);
  for (int ci = 0; ci < nc; ++ci)
    c[ci] = 0;
  for (int ni = 0; ni < nen; ++ni)
    for (int ci = 0; ci < nc; ++ci)
      c[ci] += nodeData[ni * nc + ci] * shapeValues[ni];
}

void Element::getNodeData()
{
  field->getData()->getElementData(entity,nodeData);
}

}//namespace apf

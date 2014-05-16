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

void Element::getGlobalGradients(Vector3 const& local,
                                 NewArray<Vector3>& globalGradients)
{
  Matrix3x3 J;
  parent->getJacobian(local,J);
  Matrix3x3 jinv = invert(J);
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

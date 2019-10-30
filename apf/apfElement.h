/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFELEMENT_H
#define APFELEMENT_H

#include "apfMesh.h"
#include "apfField.h"
#include "apfFieldData.h"
#include "apfShape.h"

namespace apf {

Matrix3x3 getJacobianInverse(Matrix3x3 J, int dim);

class EntityShape;
class VectorElement;

template <class T>
class ElementBase
{
public:
  ElementBase(FieldBase* f, MeshEntity* e);
  ElementBase(FieldBase* f, VectorElement* p);
  virtual ~ElementBase() {}
  void getGlobalGradients(Vector3 const& local,
                          NewArray<Vector3>& globalGradients);
  int getType() { return mesh->getType(entity); }
  int getDimension() { return Mesh::typeDimension[getType()]; }
  int getOrder() { return field->getShape()->getOrder(); }
  VectorElement* getParent() { return parent; }
  MeshEntity* getEntity() { return entity; }
  Mesh* getMesh() { return mesh; }
  EntityShape * getShape() { return shape; }
  FieldBase* getFieldBase() { return field; }
  void getComponents(Vector3 const& xi, T * c)
  {
    NewArray<double> shapeValues;
    shape->getValues(mesh, entity, xi, shapeValues);
    for (int ci = 0; ci < nc; ++ci)
      c[ci] = 0;
    for (int ni = 0; ni < nen; ++ni)
      for (int ci = 0; ci < nc; ++ci)
        c[ci] += nodeData[ni * nc + ci] * shapeValues[ni];
  }
protected:
  void init(FieldBase* f, MeshEntity* e, VectorElement* p)
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
  void getNodeData()
  {
    reinterpret_cast<FieldDataOf<T>*>(field->getData())->getElementData(entity,nodeData);
  }
  NewArray<T> nodeData;
  FieldBase* field;
  Mesh* mesh;
  MeshEntity* entity;
  EntityShape* shape;
  VectorElement* parent;
  int nen;
  int nc;
};


template <class T>
int countNodes(ElementBase<T> * e)
{
  return e->getShape()->countNodes();
}

template <class T>
void getShapeValues(ElementBase<T> * e, Vector3 const& local, NewArray<double>& values)
{
  e->getShape()->getValues(e->getMesh(),e->getEntity(),local,values);
}

template <class T>
void getShapeGrads(ElementBase<T> * e, Vector3 const& local, NewArray<Vector3>& grads)
{
  e->getGlobalGradients(local,grads);
}

}//namespace apf

#endif

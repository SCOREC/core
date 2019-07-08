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

class EntityShape;
class VectorElement;

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
protected:
  void init(FieldBase* f, MeshEntity* e, VectorElement* p);
  virtual void getNodeData() = 0;
  FieldBase* field;
  Mesh* mesh;
  MeshEntity* entity;
  EntityShape* shape;
  VectorElement* parent;
  int nen;
  int nc;
};

// not to be confused with ElementOf<T>, which is about
//  Scalar/Vector/Matrix node/value types
//  this is about the underlying scalar type (real/complex)
template <class T>
class ElementT : public ElementBase
{
public:
  ElementT(FieldBase* f, MeshEntity* e) : ElementBase(f,e) { }
  ElementT(FieldBase* f, VectorElement* p) : ElementBase(f,p) { }
  virtual ~ElementT() {}
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
  virtual void getNodeData()
  {
    reinterpret_cast<FieldDataOf<T>*>(field->getData())->getElementData(entity,nodeData);
  }
  NewArray<T> nodeData;
};

class Element : public ElementT<double>
{
  public:
    Element(FieldBase* f, MeshEntity* e)
      : ElementT<double>(f,e)
    { }
    Element(FieldBase* f, VectorElement* p)
      : ElementT<double>(f,p)
    { }
};

Matrix3x3 getJacobianInverse(Matrix3x3 J, int dim);
int countNodes(ElementBase * e);

}//namespace apf

#endif

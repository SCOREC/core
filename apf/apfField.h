/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFFIELD_H
#define APFFIELD_H

#include <string>
#include "apfMesh.h"

namespace apf {

class Element;
class FieldData;

class FieldBase
{
  public:
    void init(const char* n,
              Mesh* m,
              FieldShape* s,
              FieldData* d);
    virtual ~FieldBase();
    /* returns the number of components per node:
       1 for a scalar field and Nedelec type fields
       3 for a vector field
       9 for a matrix field, etc. */
    virtual int countComponents() const = 0;
    virtual int getScalarType() = 0;
    const char* getName() {return name.c_str();}
    Mesh* getMesh() {return mesh;}
    FieldShape* getShape() {return shape;}
    FieldData* getData() {return data;}
    int countNodesOn(MeshEntity* e);
    int countValuesOn(MeshEntity* e);
    void changeData(FieldData* d);
    void rename(const char* newName);
  protected:
    std::string name;
    Mesh* mesh;
    FieldShape* shape;
    FieldData* data;
};

template <class T>
class FieldDataOf;

class VectorElement;

class Field : public FieldBase
{
  public:
    virtual Element* getElement(VectorElement* e) = 0;
    // Think of this as the dof types, that is:
    // ScalarField will have 1 dof per node
    // VectorField will have 3 dof per noed
    // MixedVectorField (eg Nedelec) will have 1 dof per node
    virtual int getValueType() const = 0;
    // Think of this as the number of components for the shape, that is:
    // Regular shape function will have 1 value per node, but Nedelec type
    // shapes will have 3 values per node
    virtual int getShapeType() const = 0;
    virtual int getScalarType() {return Mesh::DOUBLE;}
    FieldDataOf<double>* getData();
    virtual void project(Field* from) = 0;
    virtual void axpy(double a, Field* x) = 0;
};

class FieldOp
{
  public:
    virtual bool inEntity(MeshEntity* e);
    virtual void outEntity();
    virtual void atNode(int node);
    void apply(FieldBase* f);
};

Field* makeField(
    Mesh* m,
    const char* name,
    int valueType,
    int components,
    FieldShape* shape,
    FieldData* data);

} //namespace apf

#endif

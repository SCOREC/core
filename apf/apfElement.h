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
#include "apfShape.h"

namespace apf {

class EntityShape;
class FieldShape;
class VectorElement;

class Element
{
  public:
    Element(Field* f, MeshEntity* e);
    Element(Field* f, VectorElement* p);
    virtual ~Element();
    void getGlobalGradients(Vector3 const& local,
                            NewArray<Vector3>& globalGradients);
    int getType() {return mesh->getType(entity);}
    int getDimension() {return Mesh::typeDimension[getType()];}
    int getOrder() {return field->getShape()->getOrder();}
    VectorElement* getParent() {return parent;}
    MeshEntity* getEntity() {return entity;}
    Mesh* getMesh() {return mesh;}
    EntityShape* getShape() {return shape;}
    FieldShape* getFieldShape() {return field->getShape();}
    void getComponents(Vector3 const& xi, double* c);
    void getElementNodeData(NewArray<double>& d);
  protected:
    void init(Field* f, MeshEntity* e, VectorElement* p);
    void getNodeData();
    Field* field;
    Mesh* mesh;
    MeshEntity* entity;
    EntityShape* shape;
    VectorElement* parent;
    int nen;
    int nc;
    NewArray<double> nodeData;
};

Matrix3x3 getJacobianInverse(Matrix3x3 J, int dim);

}//namespace apf

#endif

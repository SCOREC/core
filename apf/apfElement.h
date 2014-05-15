/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
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
    EntityShape* getShape() {return shape;}
    void getComponents(Vector3 const& xi, double* c);
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

}//namespace apf

#endif

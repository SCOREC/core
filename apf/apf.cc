/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apf.h"
#include "apfScalarField.h"
#include "apfScalarElement.h"
#include "apfVectorField.h"
#include "apfVectorElement.h"
#include "apfMatrixField.h"
#include "apfMatrixElement.h"
#include "apfPackedField.h"
#include "apfIntegrate.h"
#include "apfArrayData.h"
#include "apfTagData.h"
#include "apfUserData.h"
#include <cstdio>
#include <cstdlib>

namespace apf {

void destroyMesh(Mesh* m)
{
  while (m->countFields())
    destroyField(m->getField(0));
  delete m;
}

MeshEntity* getMeshEntity(MeshElement* me)
{
  return me->getEntity();
}

MeshElement* createMeshElement(Mesh* m, MeshEntity* e)
{
  return new VectorElement(static_cast<VectorField*>(
        m->getCoordinateField()), e);
}

void destroyMeshElement(MeshElement* e)
{
  delete e;
}

Field* makeField(
    Mesh* m,
    const char* name,
    int valueType,
    int components,
    FieldShape* shape,
    FieldData* data)
{
  assert( ! m->findField(name));
  Field* f = 0;
  if (valueType == SCALAR)
    f = new ScalarField();
  else if (valueType == VECTOR)
    f = new VectorField();
  else if (valueType == MATRIX)
    f = new MatrixField();
  else if (valueType == PACKED)
    f = new PackedField(components);
  else
    fail("invalid valueType in field construction\n");
  f->init(name,m,shape,data);
  m->addField(f);
  return f;
}

Field* createField(Mesh* m, const char* name, int valueType, FieldShape* shape)
{
  return makeField(m, name, valueType, 0, shape, new TagDataOf<double>);
}

Field* createLagrangeField(Mesh* m, const char* name, int valueType, int order)
{
  return createField(m,name,valueType,getLagrange(order));
}

Field* createStepField(Mesh* m, const char* name, int valueType)
{
  return createField(m,name,valueType,getConstant(m->getDimension()));
}

Field* createIPField(Mesh* m, const char* name, int valueType, int order)
{
  return createField(m,name,valueType,getIPShape(m->getDimension(),order));
}

Field* createFieldOn(Mesh* m, const char* name, int valueType)
{
  return createField(m,name,valueType,m->getShape());
}

Field* createPackedField(Mesh* m, const char* name, int components)
{
  return makeField(m, name, PACKED, components, m->getShape(),
      new TagDataOf<double>());
}

Field* cloneField(Field* f, Mesh* onto)
{
  return makeField(onto, f->getName(), f->getValueType(), f->countComponents(),
      f->getShape(), f->getData()->clone());
}

Mesh* getMesh(Field* f)
{
  return f->getMesh();
}

bool hasEntity(Field* f, MeshEntity* e)
{
  return f->getData()->hasEntity(e);
}

const char* getName(Field* f)
{
  return f->getName();
}

int getValueType(Field* f)
{
  return f->getValueType();
}

void destroyField(Field* f)
{
  if (!f)
    return;
  getMesh(f)->removeField(f);
  delete f;
}

void setScalar(Field* f, MeshEntity* e, int node, double value)
{
  ScalarField* field = static_cast<ScalarField*>(f);
  field->setNodeValue(e,node,value);
}

double getScalar(Field* f, MeshEntity* e, int node)
{
  ScalarField* field = static_cast<ScalarField*>(f);
  double value;
  field->getNodeValue(e,node,value);
  return value;
}

void setVector(Field* f, MeshEntity* e, int node, Vector3 const& value)
{
  VectorField* field = static_cast<VectorField*>(f);
  field->setNodeValue(e,node,value);
}

void getVector(Field* f, MeshEntity* e, int node, Vector3& value)
{
  VectorField* field = static_cast<VectorField*>(f);
  field->getNodeValue(e,node,value);
}

void setMatrix(Field* f, MeshEntity* e, int node, Matrix3x3 const& value)
{
  MatrixField* field = static_cast<MatrixField*>(f);
  field->setNodeValue(e,node,value);
}

void getMatrix(Field* f, MeshEntity* e, int node, Matrix3x3& value)
{
  MatrixField* field = static_cast<MatrixField*>(f);
  field->getNodeValue(e,node,value);
}

void setComponents(Field* f, MeshEntity* e, int node, double const* components)
{
  f->getData()->setNodeComponents(e,node,components);
}

void getComponents(Field* f, MeshEntity* e, int node, double* components)
{
  f->getData()->getNodeComponents(e,node,components);
}

Element* createElement(Field* f, MeshElement* e)
{
  return f->getElement(e);
}

Element* createElement(Field* f, MeshEntity* e)
{
  return new Element(f,e);
}

void destroyElement(Element* e)
{
  delete e;
}

MeshElement* getMeshElement(Element* e)
{
  return e->getParent();
}

MeshEntity* getMeshEntity(Element* e)
{
  return e->getEntity();
}

double getScalar(Element* e, Vector3 const& param)
{
  ScalarElement* element = static_cast<ScalarElement*>(e);
  return element->getValue(param);
}

void getGrad(Element* e, Vector3 const& param, Vector3& g)
{
  ScalarElement* element = static_cast<ScalarElement*>(e);
  element->grad(param,g);
}

void getVector(Element* e, Vector3 const& param, Vector3& value)
{
  VectorElement* element = static_cast<VectorElement*>(e);
  value = element->getValue(param);
}

double getDiv(Element* e, Vector3 const& param)
{
  VectorElement* element = static_cast<VectorElement*>(e);
  return element->div(param);
}

void getCurl(Element* e, Vector3 const& param, Vector3& c)
{
  VectorElement* element = static_cast<VectorElement*>(e);
  return element->curl(param,c);
}

void getVectorGrad(Element* e, Vector3 const& param, Matrix3x3& deriv)
{
  VectorElement* element = static_cast<VectorElement*>(e);
  return element->grad(param,deriv);
}

void getMatrix(Element* e, Vector3 const& param, Matrix3x3& value)
{
  MatrixElement* element = static_cast<MatrixElement*>(e);
  value = element->getValue(param);
}

void getComponents(Element* e, Vector3 const& param, double* components)
{
  e->getComponents(param,components);
}

int countIntPoints(MeshElement* e, int order)
{
  return getIntegration(e->getType())->getAccurate(order)->countPoints();
}

void getIntPoint(MeshElement* e, int order, int point, Vector3& param)
{
  IntegrationPoint const* p = 
    getIntegration(e->getType())->getAccurate(order)->getPoint(point);
  param = p->param;
}

double getIntWeight(MeshElement* e, int order, int point)
{
  IntegrationPoint const* p = 
    getIntegration(e->getType())->getAccurate(order)->getPoint(point);
  return p->weight;
}

void mapLocalToGlobal(MeshElement* e, Vector3 const& local, Vector3& global)
{
  global = e->getValue(local);
}

double getDV(MeshElement* e, Vector3 const& param)
{
  return e->getDV(param);
}

int getOrder(MeshElement* e)
{
  return e->getOrder();
}

void getJacobian(MeshElement* e, Vector3 const& local, Matrix3x3& j)
{
  e->getJacobian(local,j);
}

int countNodes(Element* e)
{
  return e->getShape()->countNodes();
}

void getScalarNodes(Element* e, NewArray<double>& values)
{
  ElementOf<double>* element = static_cast<ElementOf<double>*>(e);
  element->getValues(values);
}

void getVectorNodes(Element* e, NewArray<Vector3>& values)
{
  ElementOf<Vector3>* element = static_cast<ElementOf<Vector3>*>(e);
  element->getValues(values);
}

void getMatrixNodes(Element* e, NewArray<Matrix3x3>& values)
{
  ElementOf<Matrix3x3>* element = static_cast<ElementOf<Matrix3x3>*>(e);
  element->getValues(values);
}

void getShapeValues(Element* e, Vector3 const& local,
    NewArray<double>& values)
{
  e->getShape()->getValues(local,values);
}

void getShapeGrads(Element* e, Vector3 const& local,
    NewArray<Vector3>& grads)
{
  e->getGlobalGradients(local,grads);
}

FieldShape* getShape(Field* f)
{
  return f->getShape();
}

int countComponents(Field* f)
{
  return f->countComponents();
}

void getGaussPoint(int type, int order, int point, Vector3& param)
{
  param = getIntegration(type)->getAccurate(order)->getPoint(point)->param;
}

int countGaussPoints(int type, int order)
{
  return getIntegration(type)->getAccurate(order)->countPoints();
}

int getDimension(MeshElement* me)
{
  return me->getDimension();
}

void synchronize(Field* f, Sharing* shr)
{
  f->getData()->synchronize(shr);
}

void accumulate(Field* f, Sharing* shr)
{
  accumulateFieldData(f->getData(), shr);
}

void fail(const char* why)
{
  fprintf(stderr,"APF FAILED: %s\n",why);
  abort();
}

void freeze(Field* f)
{
  f->getMesh()->hasFrozenFields = true;
  freezeFieldData<double>(f);
}

void unfreeze(Field* f)
{
  unfreezeFieldData<double>(f);
}

bool isFrozen(Field* f)
{
  return f->getData()->isFrozen();
}

Function::~Function()
{
}

Field* createUserField(Mesh* m, const char* name, int valueType, FieldShape* s,
    Function* f)
{
  return makeField(m, name, valueType, 0, s, new UserData(f));
}

void copyData(Field* to, Field* from)
{
  copyFieldData(to->getData(), from->getData());
}

void projectField(Field* to, Field* from)
{
  to->project(from);
}

void axpy(double a, Field* x, Field* y)
{
  y->axpy(a, x);
}

}//namespace apf

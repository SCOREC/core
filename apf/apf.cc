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
#include <pcu_util.h>

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

MeshElement* createMeshElement(apf::Field* f, MeshEntity*e)
{
  PCU_DEBUG_ASSERT(apf::getValueType(f) == apf::VECTOR);
  return new VectorElement(static_cast<VectorField*>(f), e);
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
  PCU_ALWAYS_ASSERT( ! m->findField(name));
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

Field* createGeneralField(
    Mesh* m,
    const char* name,
    int valueType,
    int components,
    FieldShape* shape)
{
  return makeField(m, name, valueType, components, shape,
                   new TagDataOf<double>);
}

Field* createField(Mesh* m, const char* name, int valueType, FieldShape* shape)
{
  return createGeneralField(m, name, valueType, 0, shape);
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

Field* createPackedField(Mesh* m, const char* name, int components,
    apf::FieldShape* shape)
{
  if (!shape)
    shape = m->getShape();
  return createGeneralField(m, name, PACKED, components, shape);
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
  double tmp[1] = {value};
  field->setNodeValue(e,node,tmp);
}

double getScalar(Field* f, MeshEntity* e, int node)
{
  ScalarField* field = static_cast<ScalarField*>(f);
  double value[1];
  field->getNodeValue(e,node,value);
  return value[0];
}

void setVector(Field* f, MeshEntity* e, int node, Vector3 const& value)
{
  VectorField* field = static_cast<VectorField*>(f);
  Vector3 tmp[1] = {value};
  field->setNodeValue(e,node,tmp);
}

void getVector(Field* f, MeshEntity* e, int node, Vector3& value)
{
  VectorField* field = static_cast<VectorField*>(f);
  Vector3 tmp[1];
  field->getNodeValue(e,node,tmp);
  value = tmp[0];
}

void setMatrix(Field* f, MeshEntity* e, int node, Matrix3x3 const& value)
{
  MatrixField* field = static_cast<MatrixField*>(f);
  Matrix3x3 tmp[1] = {value};
  field->setNodeValue(e,node,tmp);
}

void getMatrix(Field* f, MeshEntity* e, int node, Matrix3x3& value)
{
  MatrixField* field = static_cast<MatrixField*>(f);
  Matrix3x3 tmp[1];
  field->getNodeValue(e,node,tmp);
  value = tmp[0];
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

void getJacobianInv(MeshElement* e, Vector3 const& local, Matrix3x3& jinv)
{
  Matrix3x3 j;
  getJacobian(e, local, j);
  jinv = getJacobianInverse(j, e->getDimension());
}

double computeCosAngle(Mesh* m, MeshEntity* pe, MeshEntity* e1, MeshEntity* e2,
    const Matrix3x3& Q)
{
  int peType = m->getType(pe);
  int e1Type = m->getType(e1);
  int e2Type = m->getType(e2);
  double cosAngle = 1.0; // assuming default value for angle is 0.0
  if (peType == Mesh::TET) {
    PCU_ALWAYS_ASSERT_VERBOSE(e1Type != Mesh::VERTEX && e2Type != Mesh::VERTEX,
    	"Cannot compute angle b/w vert and another entity. Aborting! ");

    PCU_ALWAYS_ASSERT_VERBOSE(e1Type != Mesh::TET && e2Type != Mesh::TET,
    	"e1 and e2 must be of type TRIANGLE or EDGE. Aborting! ");

    cosAngle = computeCosAngleInTet(m, pe, e1, e2, Q);
  }
  else if (peType == Mesh::TRIANGLE) {
    PCU_ALWAYS_ASSERT_VERBOSE(e1Type != Mesh::VERTEX && e2Type != Mesh::VERTEX,
    	"Cannot compute angle b/w vert and another entity. Aborting! ");

    PCU_ALWAYS_ASSERT_VERBOSE(e1Type != Mesh::TRIANGLE && e2Type != Mesh::TRIANGLE,
    	"e1 and e2 must be of type EDGE. Aborting! ");

    cosAngle = computeCosAngleInTri(m, pe, e1, e2, Q);
  } else {
    printf("The requested angle computation is not implemented. Aborting! \n");
    abort();
  }
  return cosAngle;
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
  e->getShape()->getValues(e->getMesh(), e->getEntity(), local,values);
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
  synchronizeFieldData<double>(f->getData(), shr);
}

void accumulate(Field* f, Sharing* shr, bool delete_shr)
{
  sharedReduction(f, shr, delete_shr, ReductionSum<double>() );
}

void sharedReduction(Field* f, Sharing* shr, bool delete_shr, 
            const ReductionOp<double>& reduce_op )
{
  reduceFieldData(f->getData(), shr, delete_shr, reduce_op);
}
void fail(const char* why)
{
  fprintf(stderr,"APF FAILED: %s\n",why);
  abort();
}

void freeze(Field* f)
{
  if (isFrozen(f)) return;
  f->getMesh()->hasFrozenFields = true;
  freezeFieldData<double>(f);
}

void unfreeze(Field* f)
{
  if (isFrozen(f))
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
  copyFieldData(from->getData(), to->getData());
}

void projectField(Field* to, Field* from)
{
  to->project(from);
}

void axpy(double a, Field* x, Field* y)
{
  y->axpy(a, x);
}

void renameField(Field* f, const char* name)
{
  Mesh* m = f->getMesh();
  PCU_ALWAYS_ASSERT( ! m->findField(name));
  f->rename(name);
}

/* bng: the two functions below are kind of incosistent with
   the others in this file (too long). can we stick their
   guts in a better place? */
void getBF(FieldShape* s, MeshElement* e, Vector3 const& p,
    NewArray<double>& BF)
{
  Mesh* m = e->getMesh();
  MeshEntity* ent = getMeshEntity(e);
  Mesh::Type t = m->getType(ent);
  EntityShape* es = s->getEntityShape(t);
  es->getValues(m, ent, p, BF);
}

void getGradBF(FieldShape* s, MeshElement* e, Vector3 const& p,
    NewArray<Vector3>& gradBF)
{
  Mesh* m = e->getMesh();
  MeshEntity* ent = getMeshEntity(e);
  Mesh::Type t = m->getType(ent);
  EntityShape* es = s->getEntityShape(t);
  Matrix3x3 jinv;
  getJacobianInv(e, p, jinv);
  NewArray<Vector3> gbf;
  es->getLocalGradients(m, ent, p, gbf);
  int nen = es->countNodes(); 
  gradBF.allocate(nen);
  for (int i=0; i < nen; ++i)
    gradBF[i] = jinv * gbf[i];

}

}//namespace apf

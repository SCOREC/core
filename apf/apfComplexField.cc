#include "apfComplexField.h"
#include "apfArrayData.h"
#include "apfElement.h"
#include "apfTagData.h"
#include "apfZero.h"
#include <pcu_util.h>

namespace apf
{

FieldDataOf<double_complex>* ComplexField::getData()
{
  return static_cast<FieldDataOf<double_complex>*>(data);
}

ComplexElement * ComplexPackedField::getElement(VectorElement * e)
{
  return new ComplexElement(this,e);
}

void ComplexPackedField::project(Field*)
{
  fail("ComplexPackedField::project is unimplemented");
}

void ComplexPackedField::axpy(double_complex, Field*)
{
  fail("ComplexPackedField::axpy is unimplemented");
}

ComplexField* makeComplexField(Mesh* m,
                               const char* name,
                               int valueType,
                               int components,
                               FieldShape * shape,
                               FieldData * data)
{
  PCU_ALWAYS_ASSERT(!m->findComplexField(name));
  ComplexField* f = 0;
  if(valueType == PACKED)
    f = new ComplexPackedField(components);
  else
    fail("invalid valueType in complex field construction\n");
  f->init(name,m,shape,data);
  m->addComplexField(f);
  return f;
}

ComplexField* createComplexPackedField(Mesh* m,
                                       const char* name,
                                       int valueType,
                                       int components,
                                       FieldShape* shape)
{
  return makeComplexField(m,name,valueType,components,shape,
                          new TagDataOf<double_complex>);
}

void freeze(ComplexField* f)
{
  if(isFrozen(f)) return;
  f->getMesh()->hasFrozenFields = true;
  freezeFieldData<double_complex>(f);
}

void unfreeze(ComplexField* f)
{
  if(isFrozen(f))
    unfreezeFieldData<double_complex>(f);
}

bool isFrozen(ComplexField* f)
{
  return isFrozen(static_cast<FieldBase*>(f));
}

void zero(ComplexField* f)
{
  ZeroOp<double_complex> op(f);
  op.apply(f);
}

void setComponents(ComplexField* f, MeshEntity* e, int node, double_complex const * components)
{
  PCU_DEBUG_ASSERT(f->getValueType() == SCALAR || f->getValueType() == PACKED);
  f->getData()->setNodeComponents(e,node,components);
}

void getComponents(ComplexField* f, MeshEntity* e, int node, double_complex * components)
{
  PCU_DEBUG_ASSERT(f->getValueType() == SCALAR || f->getValueType() == PACKED);
  f->getData()->getNodeComponents(e,node,components);
}

ComplexElement* createElement(ComplexField* f, MeshElement* e)
{
  return f->getElement(e);
}

ComplexElement* createElement(ComplexField* f, MeshEntity* e)
{
  return new ComplexElement(f,e);
}

void destroyElement(ComplexElement* e)
{
  delete e;
}

void getComponents(ComplexElement* e, Vector3 const& param, double_complex* components)
{
  e->getComponents(param,components);
}

int countNodes(ComplexElement* e)
{
  return countNodes(static_cast<ElementBase*>(e));
}



}

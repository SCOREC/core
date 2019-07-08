#ifndef APFCOMPLEXFIELD_H_
#define APFCOMPLEXFIELD_H_

#include "apfElement.h"
#include "apfField.h"
#include "apf.h"

namespace apf
{

class ComplexElement : public ElementT<double_complex>
{
public:
  ComplexElement(FieldBase* f, MeshEntity* e)
    : ElementT<double_complex>(f,e)
  { }
  ComplexElement(FieldBase* f, VectorElement* p)
    : ElementT<double_complex>(f,p)
  { }
};

class ComplexField : public FieldBase
{
public:
  virtual ComplexElement * getElement(VectorElement*) = 0;
  virtual int getValueType() const = 0;
  virtual int getScalarType() { return Mesh::COMPLEX; }
  FieldDataOf<double_complex>* getData();
  virtual void project(Field * frm) = 0;
  virtual void axpy(double_complex a, Field* x) = 0;
};

class ComplexPackedField : public ComplexField
{
public:
  ComplexPackedField(int c):components(c) {}
  virtual ~ComplexPackedField() {}
  virtual ComplexElement * getElement(VectorElement* e);
  virtual int getValueType() const { return PACKED; }
  virtual int countComponents() const { return components; }
  virtual void project(Field* frm);
  virtual void axpy(double_complex a, Field* x);
private:
  int components;
};

}

#endif

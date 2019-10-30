#ifndef APFCOMPLEXFIELD_H_
#define APFCOMPLEXFIELD_H_

#include "apfElement.h"
#include "apfField.h"
#include "apf.h"

namespace apf
{

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

template <class T, typename = void>
class ComplexFieldOf;

template <class T>
void project(ComplexFieldOf<T>* to, ComplexFieldOf<T>* from) {}

template <class T>
void axpy(double_complex a, ComplexFieldOf<T>* x, ComplexFieldOf<T>* y) {}

template <class T>
class ComplexFieldOf<T,typename std::enable_if<std::is_standard_layout<T>::value>::type > : public ComplexField
{
public:
  virtual ~ComplexFieldOf() {}
  void setNodeValue(MeshEntity* e, int node, T const value)
  {
    getData()->setNodeComponents(
      e,node,reinterpret_cast<double_complex const*>(value));
  }
  void getNodeValue(MeshEntity* e, int node, T* value)
  {
    getData()->getNodeComponents(
      e,node,reinterpret_cast<double_complex*>(value));
  }
  void project(ComplexField* from)
  {
    apf::project<T>(this,static_cast<ComplexFieldOf<T>*>(from));
  }
  void axpy(double_complex a, ComplexField* x)
  {
    apf::axpy<T>(this,a,reinterpret_cast<ComplexFieldOf<T>*>(x));
  }
};

template <class T>
class ComplexElementOf : public ComplexElement
{
public:
  ComplexElementOf(ComplexFieldOf<T>* f, MeshEntity* e)
    : Element(f,e)
  { }
  ComplexElementOf(ComplexFieldOf<T>* f, VectorElement* p)
    : Element(f,p)
  { }
  virtual ~ComplexElementOf() { }
  T* getNodeValues()
  {
    return reinterpret_cast<T*>(&(this->nodeData[0]));
  }
  T getValue(Vector3 const& local)
  {
    T value[1];
    getComponents(local, reinterpret_cast<double_complex*>(value));
    return value[0];
  }
  void getValues(NewArray<T>& values, int nc = 1)
  {
    values.allocate(nen * nc);
    T* nodeValues = getNodeValues();
    for(int i=0; i < nen * nc; ++i)
      values[i] = nodeValues[i];
  }
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

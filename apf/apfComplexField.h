#ifndef APFCOMPLEXFIELD_H_
#define APFCOMPLEXFIELD_H_

#include "apfField.h"
#include "apf.h"

namespace apf
{
class ComplexPackedField : public Field
{
public:
  ComplexPackedField(int c) : components(c) {}
  virtual ~ComplexPackedField() {}
  virtual Element * getElement(VectorElement*);
  virtual int getValueType() const { return COMPLEX_PACKED; }
  virtual int countComponents() const { return components; }
  virtual void project(Field  * frm);
  virtual void axpy(double a, Field * x);
private:
  int components;
};

}

#endif

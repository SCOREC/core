#include "apfComplexField.h"
#include "apfElement.h"

namespace apf {

Element * ComplexPackedField::getElement(VectorElement * e)
{
  return new Element(this,e);
}

void ComplexPackedField::project(Field*)
{
  fail("ComplexPackedField::project unimplemented");
}

void ComplexPackedField::axpy(double, Field*)
{
  fail("ComplexPackedField::axpy unimplemented");
}

}

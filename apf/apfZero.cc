#include "apfZero.h"
#include "apfComplex.h"
#include "apfComplexField.h"
namespace apf
{
template <>
void setComponents<double>(FieldBase* f, MeshEntity* e, int node, double const * components)
{
  setComponents(static_cast<Field*>(f),e,node,components);
}

template <>
void setComponents<double_complex>(FieldBase* f, MeshEntity* e, int node, double_complex const * components)
{
  setComponents(static_cast<ComplexField*>(f),e,node,components);
}

}

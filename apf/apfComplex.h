#ifndef APFCOMPLEX_H_
#define APFCOMPLEX_H_

#ifdef C_COMPLEX
  #include <complex.h>
#endif

#define CXX_COMPLEX 1
#ifdef CXX_COMPLEX
  #include <complex>
  using double_complex = std::complex<double>;
#endif

namespace apf
{

// forward decls for the interface
class ComplexElement;
class ComplexField;
class Mesh;
class FieldShape;
class MeshEntity;
class VectorElement;
typedef VectorElement MeshElement;
class Vector3;
template <class T>
class NewArray;

ComplexField* createComplexField(Mesh* m,
                                 const char* name,
                                 int valueType,
                                 int components,
                                 FieldShape* shape);

void freeze(ComplexField* f);
void unfreeze(ComplexField* f);
bool isFrozen(ComplexField* f);

/** \brief Return the contiguous array storing this field.
  \details This function is only defined for fields
  which are using array storage, for which apf::isFrozen
  returns true.
  \note If the underlying field data type is NOT double_complex,
  this will cause an assert fail in all compile modes.
 */
double_complex* getComplexArrayData(ComplexField * f);
void zeroField(ComplexField* f);

void setComponents(ComplexField* f, MeshEntity* e, int node, double_complex const * components);
void getComponents(ComplexField* f, MeshEntity* e, int node, double_complex * components);

ComplexElement* createElement(ComplexField* f, MeshElement* e);
ComplexElement* createElement(ComplexField* f, MeshEntity* e);
void destroyElement(ComplexElement* e);

MeshElement* getMeshElement(ComplexElement* e);
MeshEntity* getMeshEntity(ComplexElement* e);

void getComponents(ComplexElement* e, Vector3 const& param, double_complex* components);
int countNodes(ComplexElement* e);
void getShapeValues(ComplexElement* e);
void getShapeGrads(ComplexElement* e);
void getPackedNodes(ComplexElement* e, NewArray<double_complex>& values);

}

#endif

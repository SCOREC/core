#ifndef APFCOMPLEX_H_
#define APFCOMPLEX_H_

#include "apfComplexType.h"

namespace apf
{

// forward decls for the interface
class Mesh;
class FieldShape;
class MeshEntity;
class VectorElement;
typedef VectorElement MeshElement;
class Vector3;
template <class T>
class NewArray;

ComplexField* createComplexPackedField(Mesh* m,
                                       const char* name,
                                       int components,
                                       FieldShape* shape = NULL);

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
void getShapeValues(ComplexElement* e, Vector3 const& local, NewArray<double>& values);
void getShapeGrads(ComplexElement* e, Vector3 const& local, NewArray<double>& grades);
void getPackedNodes(ComplexElement* e, NewArray<double_complex>& values);

}

#endif

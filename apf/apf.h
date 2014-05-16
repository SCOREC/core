/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_H
#define APF_H

#include "apfMatrix.h"
#include "apfNew.h"

/** \file apf.h
  * \brief The Attached Parallel Fields interface.
  *
  * \details This is the API available to users of APF.
  */

/** \namespace apf
  * \brief All APF symbols are contained in this namespace.
  *
  * \details Wrapping a namespace over everything gives reasonable
  * insurance against future symbol conflicts with other packages
  * while very common names are being used, like Mesh or VECTOR.
  * Users will be able to use the simple names directly with
  * a using statement, likewise all APF code is written without apf::
  */
namespace apf {

class Field;
class Element;
class Mesh;
class MeshEntity;
class VectorElement;
typedef VectorElement MeshElement;
class FieldShape;

/** \brief Destroys an apf::Mesh.
  *
  * \details This only destroys the apf::Mesh object, the underlying
  * mesh database is unaffected. Mesh objects are wrappers over mesh
  * databases, and are created by functions provided outside the APF core.
  */
void destroyMesh(Mesh* m);

/** \brief Creates a Mesh Element over an entity.
  * 
  * \details A Mesh Element allows queries to the coordinate field,
  * including mapping, differential and total volume, as well as
  * gauss integration point data. A Mesh Element is also required
  * to build a Field Element.
  */
MeshElement* createMeshElement(Mesh* m, MeshEntity* e);

/** \brief Retrieve the mesh entity associated with an apf::MeshElement.
  */
MeshEntity * getMeshEntity(MeshElement * me);

/** \brief Destroys a Mesh Element.
  * 
  * \details This only destroys the apf::MeshElement object,
  * the underlying mesh entity and field data are unaffected.
  */
void destroyMeshElement(MeshElement* e);

/** \brief The type of value the field stores.
  *
  * \details The near future may bring more complex tensors.
  */
enum ValueType {
 /** \brief a single scalar value. */
  SCALAR,
 /** \brief a 3D vector value */
  VECTOR,
 /** \brief a 3x3 matrix */
  MATRIX,
 /** \brief a user-defined set of components */
  PACKED,
  VALUE_TYPES
};

/** \brief Create an apf::Field using a Lagrange distribution.
  *
  * \param m the mesh over which the field is defined
  * \param name a unique name for this field
  * \param value the type of field data
  * \param order the polynomial order of the shape functions (so far 1 or 2)
  */
Field* createLagrangeField(Mesh* m, const char* name, int valueType, int order);

/** \brief Create an apf::Field using a step distribution.
  *
  * \details A step-wise distribution is a C-1 continuous field
  * defined by one node at each element, with the field value being
  * constant over the element and discontinuous at element boundaries
  *
  * \param m the mesh over which the field is defined
  * \param name a unique name for this field
  * \param value the type of field data
  */
Field* createStepField(Mesh* m, const char* name, int valueType);

/** \brief Create an apf::Field of integration point data.
  *
  * \param m the mesh over which the field is defined
  * \param name a unique name for this field
  * \param value the type of field data
  * \param order the polynomial order of accuracy to determine the integration points.
  */
Field* createIPField(Mesh* m, const char* name, int valueType, int order);

/** \brief Create an apf::Field from any built-in or user-define apf::FieldShape.
  */
Field* createField(Mesh* m, const char* name, int valueType, FieldShape* shape);

Field* createFieldOn(Mesh* m, const char* name, int valueType);

Field* createPackedField(Mesh* m, const char* name, int components);

/** \brief Retrieve the Mesh over which a Field is defined.
  */
Mesh* getMesh(Field* f);

/** \brief Returns true iff an entity has data associate with a field.
  */
bool hasEntity(Field* f, MeshEntity* e);

/** \brief Get the name of a Field.
  *
  * \details Both for use convenience and for technical reasons
  * related to tagging, each Field should have a unique name.
  */
const char* getName(Field* f);

/** \brief Retrieve the type of value a field distributes
  */
int getValueType(Field* f);

/** \brief Destroy an apf::Field.
  *
  * \details This function will also remove any field data that
  * this field attached to its Mesh domain.
  */
void destroyField(Field* f);

/** \brief Set a nodal value of a scalar field.
  *
  * \param node The node number within the entity. So far, it is just
  * 0 for vertices and for edges in 2nd order. higher order
  * will bring more nodes per edge and onto faces and such.
  */
void setScalar(Field* f, MeshEntity* e, int node, double value);

/** \brief Get the node value of a scalar field.
  *
  * \returns The field value at that node.
  */
double getScalar(Field* f, MeshEntity* e, int node);

/** \brief Set the nodal value of a vector field.
  *
  * \param value The vector value.
  */
void setVector(Field* f, MeshEntity* e, int node, Vector3 const& value);

/** \brief Get the nodal value of a vector field.
  */
void getVector(Field* f, MeshEntity* e, int node, Vector3& value);

/** \brief Set the nodal value of a matrix field.
  *
  * \param value The vector value.
  */
void setMatrix(Field* f, MeshEntity* e, int node, Matrix3x3 const& value);

/** \brief Get the nodal value of a matrix field.
  */
void getMatrix(Field* f, MeshEntity* e, int node, Matrix3x3& value);

void setComponents(Field* f, MeshEntity* e, int node, double const* components);
void getComponents(Field* f, MeshEntity* e, int node, double* components);

/** \brief Create a Field Element from a Mesh Element.
  *
  * \details A Field Element object caches elemental data for
  * use in evaluation, mapping, and integration.
  * use destroyElement to free this data.
  *
  * \param f The field which the Element will represent
  * \returns The new field Element
  */
Element* createElement(Field* f, MeshElement* e);

Element* createElement(Field* f, MeshEntity* e);

/** \brief Destroy a Field Element.
 */
void destroyElement(Element* e);

/** \brief Get the Mesh Element of a Field Element.
  *
  * \details Each Field Element operates over
  * a Mesh Element, which must be maintained
  * as long as the Field Element exists. Multiple
  * Field Elements may share a Mesh Element.
  */
MeshElement* getMeshElement(Element* e);
MeshEntity* getMeshEntity(Element* e);

/** \brief Evaluate a scalar field at a point.
  *
  * \param param The local coordinates in the element.
  * \returns The field value at that point.
  */
double getScalar(Element* e, Vector3 const& param);

/** \brief Evaluate the gradient of a scalar field at a point.
  *
  * \param param The local coordinates in the element.
  * \param grad The gradient vector at that point.
  */
void getGrad(Element* e, Vector3 const& param, Vector3& grad);

/** \brief Evaluate a vector field at a point.
  *
  * \param param The local coordinates in the element.
  * \param value The field value at that point.
  */
void getVector(Element* e, Vector3 const& param, Vector3& value);

/** \brief Evaluate the divergence of a vector field at a point.
  *
  * \param param The local coordinates in the element.
  * \returns The divergence at that point.
  */
double getDiv(Element* e, Vector3 const& param);

/** \brief Evaluate the curl of a vector field at a point.
  *
  * \param param The local coordinates in the element.
  * \param curl The curl vector at that point.
  */
void getCurl(Element* e, Vector3 const& param, Vector3& curl);

/** \brief Evaluate the gradient of a vector field at a point.
  *
  * \param param The local coordinates in the element.
  * \param deriv The gradient matrix at that point.
  */
void getVectorGrad(Element* e, Vector3 const& param, Matrix3x3& deriv);

/** \brief Evaluate the value of a matrix field.
  *
  * \param param The local coordinates in the element.
  * \param value The field value at that point.
  */
void getMatrix(Element* e, Vector3 const& param, Matrix3x3& value);

void getComponents(Element* e, Vector3 const& param, double* components);

/** \brief Get the number of integration points for an element.
  *
  * \param order the polynomial order of accuracy desired for the integration
  * (not to be confused with the polynomial order of shape functions).
  */
int countIntPoints(MeshElement* e, int order);

/** \brief Get an integration point in an element.
  *
  * \param order The polynomial order of accuracy.
  * \param point The integration point number.
  * \param param The resulting local coordinate of the integration point.
  */
void getIntPoint(MeshElement* e, int order, int point, Vector3& param);

/** \brief Get the weight of an integration point in an element.
  *
  * \details All integration point tables in APF are scaled
  * such that the sum of the weights equals the area of
  * of the parent element.
  */
double getIntWeight(MeshElement* e, int order, int point);

/** \brief Map a local coordinate to a global coordinate.
  */
void mapLocalToGlobal(MeshElement* e, Vector3 const& local, Vector3& global);

/** \brief Get the differential volume at a point.
  * 
  * \details This function is meant to provide the differential
  * measure of an element at a point, based on
  * the Jacobian determinant in the case of regions, and equivalent
  * terms for lower dimensions.
  *
  * \returns The differential volume
  */
double getDV(MeshElement* e, Vector3 const& param);

/** \brief A virtual base for user-defined integrators.
  *
  * \details Users of APF can define an Integrator object to handle
  * integrating expressions over elements and over meshes.
  * Users specify the accuracy of the integration and provide
  * accumulation callbacks which APF uses at each integration
  * point. The APF-provided process functions will perform the
  * integration over an element or mesh using the callbacks.
  * In parallel, users must provide a reduction callback
  * to turn locally accumulated values into a globally integrated
  * value.
  */
class Integrator
{
  public:
    /** \brief Construct an Integrator given an order of accuracy. */
    Integrator(int o);
    virtual ~Integrator();
    /** \brief Run the Integrator over the local Mesh. */
    void process(Mesh* m);
    /** \brief Run the Integrator over a Mesh Element. */
    void process(MeshElement* e);
    /** \brief User callback: element entry.
      *
      * \details APF will call this function every time the
      * Integrator begins operating over a new element.
      * Users can then construct Field Elements, for example.
      */
    virtual void inElement(MeshElement*) {}
    /** \brief User callback: element exit.
      *
      * \details APF will call this function once an Integrator
      * is done operating over an element. This can be used
      * to destroy Field Elements, for example.
      */
    virtual void outElement() {}
    /** \brief User callback: accumulation.
      *
      * \details APF will call this function at each integration
      * point. Users should evaluate their expression and accumulate
      * the value.
      *
      * \param p The local coordinates of the point.
      * \param w The integration weight of the point.
      * \param dV The differential volume at that point.
      */
    virtual void atPoint(Vector3 const& p, double w, double dV) = 0;
    /** \brief User callback: parallel reduction.
      *
      * \details This function should use communication to reduce
      * process-local integrations into a global mesh integration,
      * if that is the user's goal.
      */
    virtual void parallelReduce() {}
  protected:
    int order;
};

/** \brief Measures the volume, area, or length of a Mesh Element.
  *
  * \details By integrating the differential volume over the element,
  * a general measure is obtained. This correctly measures curved
  * meshes.
  *
  * \returns The measure of the element
  */
double measure(MeshElement* e);

/** \brief Returns the polynomial order of the coordinate field.
  */
int getOrder(MeshElement* e);

/** \brief Returns the Jacobian at a local point
  */
void getJacobian(MeshElement* e, Vector3 const& local, Matrix3x3& j);

/** \brief Returns the number of element nodes.
  *
  * \details This is the number of nodes affecting an
  * element, as opposed to the nodes tagged to an entity.
  */
int countNodes(Element* e);

/** \brief Returns the element nodal values for a scalar field
  */
void getScalarNodes(Element* e, NewArray<double>& values);

/** \brief Returns the element nodal values for a vector field
  */
void getVectorNodes(Element* e, NewArray<Vector3>& values);

/** \brief Returns the element nodal values for a matrix field
  */
void getMatrixNodes(Element* e, NewArray<Matrix3x3>& values);

/** \brief Returns the shape function values at a point
  */
void getShapeValues(Element* e, Vector3 const& local,
    NewArray<double>& values);

/** \brief Returns the shape function gradients at a point
  *
  * \details these are gradients with respect to global coordinates.
  */
void getShapeGrads(Element* e, Vector3 const& local,
    NewArray<Vector3>& grads);

/** \brief Retrieve the apf::FieldShape used by a field
  */
FieldShape* getShape(Field* f);

/** \brief Count the number of scalar components in the field's value type
  */
int countComponents(Field* f);

#define APF_ITERATE(t,w,i) \
for (t::iterator (i) = (w).begin(); \
     (i) != (w).end(); ++(i))

/** \brief Write a set of parallel VTK Unstructured Mesh files from an apf::Mesh
  */
void writeVtkFiles(const char* prefix, Mesh* m);
void writeOneVtkFile(const char* prefix, Mesh* m);

void getGaussPoint(int type, int order, int point, Vector3& param);
int countGaussPoints(int type, int order);

/* a special function taking into account the various
   formulae for differential volume at each dimension. */
double getJacobianDeterminant(Matrix3x3 const& J, int dimension);

int getDimension(MeshElement* me);

void synchronize(Field* f);

void accumulate(Field* f);

void fail(const char* why);

void freeze(Field* f);

void unfreeze(Field* f);

bool isFrozen(Field* f);

double* getArrayData(Field* f);
}

#endif

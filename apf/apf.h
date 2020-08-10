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
#include "apfDynamicArray.h"

#include <vector>
#include <limits>

/** \file apf.h
  * \brief The APF Field interface
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
/** \brief Mesh Elements represent the mesh coordinate vector field. */
typedef VectorElement MeshElement;
class FieldShape;
struct Sharing;
template <class T> class ReductionOp;
template <class T> class ReductionSum;

/** \brief Base class for applying operations to make a Field consistent
  * in parallel
  * \details This function gets applied pairwise to the Field values
  * from every partition, resulting in a single unique value.  No guarantees
  * are made about the order in which this function is applied to the
  * values.
  */
template <class T>
class ReductionOp
{
  public:
    /* \brief apply operation, returning a single value */
    virtual T apply(T val1, T val2) const = 0;

    /* \brief returns a value such that apply(val, neutral_val) == val */
    virtual T getNeutralElement() const = 0;
};

template <class T>
class ReductionSum : public ReductionOp<T>
{
  T apply(T val1, T val2) const { return val1 + val2; };
  T getNeutralElement() const { return 0; };
};

template <class T>
class ReductionMin : public ReductionOp<T>
{
  T apply(T val1, T val2) const { return ( (val1 < val2) ? val1 : val2 ); };
  T getNeutralElement() const { return std::numeric_limits<T>::max(); };
};

template <class T>
class ReductionMax : public ReductionOp<T>
{
  T apply(T val1, T val2) const { return ( (val1 < val2) ? val2 : val1 ); };
  T getNeutralElement() const { return std::numeric_limits<T>::min(); };
};

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

/** \brief Creates a non-standard Mesh Element over an entity.
  *
  * \details Warning: most users should use the above function unless
  * you know this is the right one for you. This function is useful for
  * creating a mesh element over a coordinate field c that may be
  * distinct from the mesh's underlying coordinate field. This may
  * be useful for updated-Lagrangian simulations, where quantities
  * must be computed on reference and updated configurations.
  */
MeshElement* createMeshElement(Field* c, MeshEntity* e);

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
 /** \brief placeholder used to set array sizes */
  VALUE_TYPES
};

/** \brief Create an apf::Field using a Lagrange distribution.
  *
  * \param m the mesh over which the field is defined
  * \param name a unique name for this field
  * \param valueType the type of field data
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
  * \param valueType the type of field data
  */
Field* createStepField(Mesh* m, const char* name, int valueType);

/** \brief Create an apf::Field of integration point data.
  *
  * \param m the mesh over which the field is defined
  * \param name a unique name for this field
  * \param valueType the type of field data
  * \param order polynomial order of accuracy
  */
Field* createIPField(Mesh* m, const char* name, int valueType, int order);

/** \brief Create a Field from any builtin or user defined FieldShape.
  */
Field* createField(Mesh* m, const char* name, int valueType, FieldShape* shape);

/** \brief Create a field using the mesh's coordinate nodal distribution */
Field* createFieldOn(Mesh* m, const char* name, int valueType);

/** \brief Create a field of N components without a tensor type.
  \details Packed fields are used to interface with applications
  whose fields are not easily categorized as tensors of order 0,1,2.
  They contain enough information to interpolate values in an
  element and such, but some higher-level functionality is left out.
  */
Field* createPackedField(Mesh* m, const char* name, int components,
    apf::FieldShape* shape = 0);

/** \brief Encompasses both packed and typed fields */
Field* createGeneralField(
    Mesh* m,
    const char* name,
    int valueType,
    int components,
    FieldShape* shape);

/** \brief Declare a copy of a field on another apf::Mesh
   \details This will just make a Field object with the same
   properties, but not fill in any data. */
Field* cloneField(Field* f, Mesh* onto);

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

/** \brief Set the nodal value from an array of component values. */
void setComponents(Field* f, MeshEntity* e, int node, double const* components);

/** \brief Copy the nodal value into an array of component values.
  \details the output array must already be allocated, apf::countComponents
  can help with this.
 */
void getComponents(Field* f, MeshEntity* e, int node, double* components);

/** \brief Create a Field Element from a Mesh Element.
  *
  * \details A Field Element object caches elemental data for
  * use in evaluation, mapping, and integration.
  * use destroyElement to free this data.
  *
  * \param f The field which the Element will represent
  * \param e An existing MeshElement for the desired entity
  * \returns The new field Element
  */
Element* createElement(Field* f, MeshElement* e);

/** \brief Create a Field Element without a parent Mesh Element
    \details Warning: most users should call the version
    which takes a MeshElement as input. Only call this
    function if you know the other one isn't right for you. */
Element* createElement(Field* f, MeshEntity* e);

/** \brief Destroy a Field Element.
 */
void destroyElement(Element* e);

/** \brief Get the Mesh Element of a Field Element.
  *
  * \details Each apf::Element operates over
  * an apf::MeshElement, which must be maintained
  * as long as the apf::Element exists. Multiple
  * apf::Elements may share an apf::MeshElement.
  */
MeshElement* getMeshElement(Element* e);

/** \brief Retrieve the mesh entity of an apf::Element. */
MeshEntity* getMeshEntity(Element* e);

/** \brief Evaluate a scalar field at a point.
  *
  * \param param The local coordinates in the element.
  * \returns The field value at that point.
  */
double getScalar(Element* e, Vector3 const& param);

/** \brief Get the gradient of a scalar field w.r.t. global coordinates
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

/** \brief Get the gradient of a vector field w.r.t. global coordinates
  *
  * \param param The local coordinates in the element.
  * \param deriv The gradient matrix at that point.
  * \details Note: The return parameter component deriv[j][i] stores the
  * value $(grad u)_{i,j} = \frac{\partial u_i} / \frac{\partial x_j}$,
  * where $u$ is the vector field of interest, i is the vector index,
  * and j is the derivative index.
  */
void getVectorGrad(Element* e, Vector3 const& param, Matrix3x3& deriv);

/** \brief Evaluate the value of a matrix field.
  *
  * \param param The local coordinates in the element.
  * \param value The field value at that point.
  */
void getMatrix(Element* e, Vector3 const& param, Matrix3x3& value);

/** \brief get the gradient of a matrix field
 * \param param The local coordinate in the element.
 * \param value the gradient of the field at a point. If the matrix is defined by A_{ij},
 * the gradient $\frac{\partial A_{ij}}{\partial X_k}=value[i*3+j+9*k]$ .
 */
void getMatrixGrad(Element* e, Vector3 const& param, Vector<27>& value);

/** \brief Evaluate a field into an array of component values. */
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
    /** \brief Run the Integrator over the local Mesh.
     * \param m mesh to integrate over
     * \param dim optional dimension to integrate over. This defaults to
     * integration over the mesh dimesion which may not be correct e.g. in the case
     * of a 1D element embeded in 3D space.
     * */
    void process(Mesh* m, int dim=-1);
    /** \brief Run the Integrator over a Mesh Element. */
    void process(MeshElement* e);
    /** \brief User callback: element entry.
      *
      * \details APF will call this function every time the
      * Integrator begins operating over a new element.
      * Users can then construct Field Elements, for example.
      */
    virtual void inElement(MeshElement*);
    /** \brief User callback: element exit.
      *
      * \details APF will call this function once an Integrator
      * is done operating over an element. This can be used
      * to destroy Field Elements, for example.
      */
    virtual void outElement();
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
    virtual void parallelReduce();
  protected:
    int order;
    int ipnode;
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

/** \brief Measures the volume, area, or length of a Mesh Entity.
  *
  * \returns The measure of the element
  */
double measure(Mesh* m, MeshEntity* e);

/** \brief Returns the polynomial order of the coordinate field.
  */
int getOrder(MeshElement* e);

/** \brief Returns the Jacobian at a local point
  */
void getJacobian(MeshElement* e, Vector3 const& local, Matrix3x3& j);

/** \brief Returns the Jacobian inverse at a local point
  *
  * \details computes the Moore-Penrose psuedo-inverse of J. This is needed
  * to handle non-square Jacobians that arise from edges embedded in 2D or
  * higher and faces embeded in 3D. If J is invertible, the Moore-Penrose
  * inverse equals the inverse.
  */
void getJacobianInv(MeshElement* e, Vector3 const& local, Matrix3x3& jinv);

/** \brief Returns the cosine of the angle between 2 entities of
  * the parent entity.
  *
  * \param m	mesh
  * \param pe	parent entity
  * \param e1	first entity
  * \param e2	second entity
  * \param Q	metric
  *
  * \details All the angles are computed based on the Jacobian, and
  * therefore, this function will work for curved elements as well. For
  * anisotropic meshes the user can provide Q (he metric value of the parent)
  * to get the angles in the metric space. Otherwise and Identity matrix must
  * be given for Q.
  */
double computeCosAngle(Mesh* m, MeshEntity* pe, MeshEntity* e1, MeshEntity* e2,
    const Matrix3x3& Q);

double computeCosAngleInTri(Mesh* m, MeshEntity* tri,
    MeshEntity* e1, MeshEntity* e2, const Matrix3x3& Q);

double computeCosAngleInTet(Mesh* m, MeshEntity* tet,
    MeshEntity* e1, MeshEntity* e2, const Matrix3x3& Q);

Vector3 computeFaceNormalAtEdgeInTet(Mesh* m, MeshEntity* tet,
    MeshEntity* face, MeshEntity* edge, Matrix3x3 Q);

Vector3 computeFaceNormalAtVertex(Mesh* m, MeshEntity* face,
    MeshEntity* vert, const Matrix3x3& Q);

Vector3 computeEdgeTangentAtVertex(Mesh* m, MeshEntity* edge,
    MeshEntity* vert, const Matrix3x3& Q);

/** \brief Returns the shortest height in a tet
  *
  * \param m	mesh
  * \param tet  tet
  * \param Q	metric (default Identity)
  *
  */
double computeShortestHeightInTet(Mesh* m, MeshEntity* tet,
    const Matrix3x3& Q = Matrix3x3(1.,0.,0.,0.,1.,0.,0.,0.,1.));

/** \brief Returns the largest height in a tet
  *
  * \param m	mesh
  * \param tet  tet
  * \param Q	metric (default Identity)
  *
  */
double computeLargestHeightInTet(Mesh* m, MeshEntity* tet,
    const Matrix3x3& Q = Matrix3x3(1.,0.,0.,0.,1.,0.,0.,0.,1.));


/** \brief Returns the shortest height in a tri
  *
  * \param m	mesh
  * \param tri  tri
  * \param Q	metric (default Identity)
  *
  */
double computeShortestHeightInTri(Mesh* m, MeshEntity* tri,
    const Matrix3x3& Q = Matrix3x3(1.,0.,0.,0.,1.,0.,0.,0.,1.));

/** \brief Returns the largest height in a tri
  *
  * \param m	mesh
  * \param tri  tri
  * \param Q	metric (default Identity)
  *
  */
double computeLargestHeightInTri(Mesh* m, MeshEntity* tri,
    const Matrix3x3& Q = Matrix3x3(1.,0.,0.,0.,1.,0.,0.,0.,1.));


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

/** \brief Returns the vector shape function values at a point
 *  \details used only for Nedelec shapes
 *  (Piola transformation used to map from parent to physical coordinates)
  */
void getVectorShapeValues(Element* e, Vector3 const& local,
    NewArray<Vector3>& values);

/** \brief Returns the vector curl shape function values at a point
 *  \details used only for Nedelec shapes
 *  (Piola transformation used to map from parent to physical coordinates)
  */
void getCurlShapeValues(Element* e, Vector3 const& local,
    NewArray<Vector3>& values);

/** \brief Retrieve the apf::FieldShape used by a field
  */
FieldShape* getShape(Field* f);

/** \brief Count the number of scalar components in the field's value type
  */
int countComponents(Field* f);

/** \brief Stress-free iteration over STL-like containers.
  \param t The type of the container
  \param w The container itself
  \param i The name to give to the iterator
  \details given an STL container or anything that implements
  Type::iterator, Type::begin, Type::end, and Type::iterator::operator++(),
  this macro fills in the boilerplate of a for() loop over this container.
 */
#define APF_ITERATE(t,w,i) \
for (t::iterator i = (w).begin(); \
     (i) != (w).end(); ++(i))

/** \brief APF_ITERATE for const containers. */
#define APF_CONST_ITERATE(t,w,i) \
for (t::const_iterator i = (w).begin(); \
     (i) != (w).end(); ++(i))

/** \brief Write a set of parallel VTK Unstructured Mesh files from an apf::Mesh
  * with binary (base64) encoding and zlib compression (if LION_COMPRESS=ON)
  * \details Nodal fields whose shape differs from the mesh shape will
  * not be output. Fields with incomplete data will not be output.
  */
void writeVtkFiles(const char* prefix, Mesh* m, int cellDim = -1);

/** \brief Write a set of parallel VTK Unstructured Mesh files from an apf::Mesh
  * with binary (base64) encoding and zlib compression (if LION_COMPRESS=ON)
  * \details Only fields whose name appears in the vector writeFields will be
  * output. Nodal fields whose shape differs from the mesh shape will not be
  * output. Fields with incomplete data will not be output.
  */
void writeVtkFiles(const char* prefix, Mesh* m,
    std::vector<std::string> writeFields, int cellDim = -1);

/** \brief Output just the .vtu file with ASCII encoding for this part.
  \details this function is useful for debugging large parallel meshes.
  */
void writeOneVtkFile(const char* prefix, Mesh* m);

/** \brief Write a set of parallel VTK Unstructured Mesh files from an apf::Mesh
  * with ASCII encoding
  * \details Nodal fields whose shape differs from the mesh shape will
  * not be output. Fields with incomplete data will not be output.
  */
void writeASCIIVtkFiles(const char* prefix, Mesh* m);

/** \brief Write a set of parallel VTK Unstructured Mesh files from an apf::Mesh
  * with ASCII encoding
  * \details Only fields whose name appears in the vector writeFields will be
  * output. Nodal fields whose shape differs from the mesh shape will not be
  * output. Fields with incomplete data will not be output.
  */
void writeASCIIVtkFiles(const char* prefix, Mesh* m,
    std::vector<std::string> writeFields);

/** \brief Output .vtk files with ASCII encoding for this part.
  \details this function is useful for debugging meshes with Nedelec 
  fields on them.
  */
void writeNedelecVtkFiles(const char* prefix, Mesh* m);

/** \brief Return the location of a gaussian integration point.
  \param type the element type, from apf::Mesh::getType
  \param order the order of the integration rule
  \param point which point of the integration rule
  \param param the resulting parent element coordinates
 */
void getGaussPoint(int type, int order, int point, Vector3& param);

/** \brief Return the number of Gaussian integration points. */
int countGaussPoints(int type, int order);

/** \brief Return the Jacobian determinant or dimensional equivalent.
  \details this is a special function taking into account the various
  formulae for differential volume, area, lenght, etc.
  \param J Jacobian matrix, vector gradient of coordinate field
           with respect to parent element coordinates
  \param dimension spacial dimension of the entity being measured
  */
double getJacobianDeterminant(Matrix3x3 const& J, int dimension);

/** \brief Return the dimension of a MeshElement's MeshEntity. */
int getDimension(MeshElement* me);

/** \brief Synchronize field values along partition boundary.
  \details Using the ownership and copies described by an apf::Sharing
  object, copy values from the owned nodes to their copies,
  possibly assigning them values for the first time.
  */
void synchronize(Field* f, Sharing* shr = 0);

/** \brief Add field values along partition boundary.
  \details Using the copies described by
  an apf::Sharing object, add up the field values of
  all copies of an entity and assign the sum as the
  value for all copies.
  */
void accumulate(Field* f, Sharing* shr = 0, bool delete_shr = false);

/** \brief Apply a reudction operator along partition boundaries
  \details Using the copies described by an apf::Sharing object, applied
  the specified operation pairwise to the values of the field on each
  partition.  No guarantee is made about hte order of the pairwise
  application
  */
void sharedReduction(Field* f, Sharing* shr, bool delete_shr,
           const ReductionOp<double>& sum = ReductionSum<double>());

/** \brief Checks whether a Field/Numbering/GlobalNumbering is complete and
 * therefore printable to visualization files.  This is a collective operation.
 */
bool isPrintable(Field* f);

/** \brief Declare failure of code inside APF.
  \details This function prints the string as an APF
  failure to stderr and then calls abort.
  It can be called from code that is part of the
  apf namespace, but not outside of that.
  */
void fail(const char* why) __attribute__((noreturn));

/** \brief Convert a Field from Tag to array storage. */
void freeze(Field* f);

/** \brief Convert a Field from array to Tag storage. */
void unfreeze(Field* f);

/** \brief Returns true iff the Field uses array storage. */
bool isFrozen(Field* f);

/** \brief Return the contiguous array storing this field.
  \details This function is only defined for fields
  which are using array storage, for which apf::isFrozen
  returns true.
 */
double* getArrayData(Field* f);

/** \brief Initialize all nodal values with all-zero components */
void zeroField(Field* f);

/** \brief User-defined Analytic Function. */
struct Function
{
  /** \brief Possible user-defined cleanup */
  virtual ~Function();
  /** \brief The user's analytic function call.
    \details For simplicity, this
    currently only supports one node per entity.
    \param e the entity on which the node is
    \param result the field component values for that node
    */
  virtual void eval(MeshEntity* e, double* result) = 0;
};

/* \brief Create a Field from a user's analytic function.
 * \details This field will use no memory, has values on all
 * nodes, and calls the user Function for all value queries.
 * Writing to this field does nothing.
 * \warning if you copy the mesh with apf::convert you may want to use
 * apf::updateUserField to update function that this field refers to. This is
 * extremely important if the analytic function you use references user fields.
 */
Field* createUserField(Mesh* m, const char* name, int valueType, FieldShape* s,
    Function* f);

/* \brief update the analytic function from a user field
 * \details this field updates the apf::Function which the UserField refers to
 */
void updateUserField(Field* field, Function* newFunc);

/** \brief Compute a nodal gradient field from a nodal input field
  \details given a nodal field, compute approximate nodal gradient
  values by giving each node a volume-weighted average of the
  gradients computed at each element around it. */
Field* recoverGradientByVolume(Field* f);

void copyData(Field* to, Field* from);

/** \brief Project a field from an existing field */
void projectField(Field* to, Field* from);

void axpy(double a, Field* x, Field* y);

void renameField(Field* f, const char* name);

/** \brief Get the basis functions over a mesh element
  *
  * \param s the field shape that defines the basis functions
  * \param e the mesh element over which the basis functions are defined
  * \param p the local coordinates at which the basis functions are evaluated
  * \param BF nodal array of basis functions evaluated at p
  */
void getBF(FieldShape* s, MeshElement* e, Vector3 const& p,
    NewArray<double>& BF);

/** \brief Get global gradients of basis functions over a mesh element
  *
  * \param gradBF nodal array of gradients of basis functions evaluated at p
  */
void getGradBF(FieldShape* s, MeshElement* e, Vector3 const& p,
    NewArray<Vector3>& gradBF);

}

#endif

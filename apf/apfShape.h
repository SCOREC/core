/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFSHAPE_H
#define APFSHAPE_H

/** \file apfShape.h
  \brief Field node distribution and shape functions */

#include "apf.h"
#include "apfNew.h"
#include <sstream>

namespace apf {

class Mesh;
class MeshEntity;

/** \brief Shape functions over this element */
class EntityShape
{
  public:
    virtual ~EntityShape();
/** \brief evaluate element shape functions
 \param xi the parent element coordinates
 \param values each entry is the shape function value for one node */
    virtual void getValues(
        Mesh* m,
        MeshEntity* e,
        Vector3 const& xi,
        NewArray<double>& values) const = 0;
/** \brief evaluate element shape function gradients
 \param xi parent element coordinates
 \param grads each entry is the shape function gradient with
              respect to parent element coordinates for one node */
    virtual void getLocalGradients(
        Mesh* m,
        MeshEntity* e,
        Vector3 const& xi,
        NewArray<Vector3>& grads) const = 0;
/** \brief evaluate element vector shape functions
 \details this is used only for Nedelec
 \param xi the parent element coordinates
 \param values each entry is the vector shape function value for one node */
    virtual void getVectorValues(
        Mesh* m,
        MeshEntity* e,
        Vector3 const& xi,
        NewArray<apf::Vector3>& values) const = 0;
/** \brief evaluate element vector curl shape functions
 \details this is used only for Nedelec
 \param xi the parent element coordinates
 \param values each entry is the vector shape function value for one node */
    virtual void getLocalVectorCurls(
        Mesh* m,
        MeshEntity* e,
        Vector3 const& xi,
        NewArray<apf::Vector3>& values) const;
/** \brief return the number of nodes affecting this element
    \details in a linear mesh, there are two nodes affecting
             and edge, three nodes affecting a triangle,
             four for a tet, etc. */
    virtual int countNodes() const = 0;
/** \brief convert from shared node order to local element order
    \details
    when two elements share nodes on a common entity,
    that entity will have its own unique orientation.
    when there is more than one node on that entity,
    this function will be queried to determine in what order
    the shared nodes go in the list of all element nodes.

    the order array should contain numbers 0 through n-1,
    where n is the number of nodes in the shared entity only.
    the shared nodes will be reordered such that node i
    moves to position order[i], and then they will be
    added to the list of element nodes as a group.

    users are encouraged to implement this with the
    help of apf::getAlignment */
    virtual void alignSharedNodes(Mesh* m,
        MeshEntity* elem, MeshEntity* shared, int order[]);
};

/** \brief Describes field distribution and shape functions
  \details these classes are typically singletons, one for
  each shape function scheme */
class FieldShape
{
  public:
    virtual ~FieldShape();
/** \brief Get the sub-descriptor for this entity type
  \param type select from apf::Mesh::Type */
    virtual EntityShape* getEntityShape(int type) = 0;
/** \brief Return true iff there are nodes on entities of this dimension
  \details this is used to skip dimensions in loops */
    virtual bool hasNodesIn(int dimension) = 0;
/** \brief Return the number of nodes associated with an entity of
           this type
    \details in a linear mesh, nodes are associated with vertices
             but there are no nodes associated with other entities.
    \param type select from apf::Mesh::Type */
    virtual int countNodesOn(int type) = 0;
/** \brief Return the polynomial order of the shape functions
  \details this is not always applicable */
    virtual int getOrder() = 0;
/** \brief Get the parent element coordinates of an element node
  \param type element type, select from apf::Mesh::Type
  \param node index from element node ordering
  \param xi parent element coordinates */
    virtual void getNodeXi(int type, int node, Vector3& xi);
/** \brief Get the parent element tangent vector of an element node
  \param type element type, select from apf::Mesh::Type
  \param node index from element node ordering
  \param t parent element tangent
  \details this is only applicable for vector shapes */
    virtual void getNodeTangent(int type, int node, Vector3& t);
/** \brief Returns true if the shape functions are vectors */
    virtual bool isVectorShape();
/** \brief Get a unique string for this shape function scheme */
    virtual const char* getName() const = 0;
    void registerSelf(const char* name);
};

/** \brief Get the Lagrangian shape function of some polynomial order
 \details we have only first and second order so far */
FieldShape* getLagrange(int order);

/** \brief Get the Serendipity shape functions of second order */
FieldShape* getSerendipity();

/** \brief Get the Constant shape function over some dimension
 \details this pseudo-shape function places a node on every element
          of the given dimension. Dimensions up to 3 are available */
FieldShape* getConstant(int dimension);
/** \brief Get the Integration Point distribution
  \param dimension The dimensionality of the elements
  \param order The order of accuracy, determines the integration points
  \details 
  This allows users to create a field which has values at the integration
  points of elements.
  Orders 1 to 3 for dimension 2 or 3 are available */
FieldShape* getIPShape(int dimension, int order);
/** \brief Get the Voronoi shape function
  \details the Voronoi FieldShape is equivalent to the IPShape except
           that it is capable of evaluating as a shape function whose
           value at any point in the element is the value of the closest
           integration point in that element. */
FieldShape* getVoronoiShape(int dimension, int order);
/** \brief Get the IP Fit shape function
  \details the IP Fit FieldShape is equivalent to the IPShape except
           that it is capable of evaluating as a shape function whose
           value at any point in the element is a polynomial fit to
           the integration point data in that element. */
FieldShape* getIPFitShape(int dimension, int order);

/** \brief Get the quadratic hierarchic shape function
  \details only first and second order so far
 */
FieldShape* getHierarchic(int order);

/** \brief Get the Nedelec shape function of a given order
  \details TODO: complete later
 */
FieldShape* getNedelec(int order);

/** \brief Get the L2 shapes of a given order and entity type
  \details
 */
FieldShape* getL2Shape(int order, int type);

/** \brief Get the H1 shapes of a given order
  \details These are hierarchic shapes that are compatible with MFEM's impl.
 */
FieldShape* getH1Shape(int order);

/** \brief Project a hierarchic field */
void projectHierarchicField(Field* to, Field* from);

/** \brief Project a Nedelec field */
void projectNedelecField(Field* to, Field* from);

/** \brief Project a L2 field */
void projectL2Field(Field* to, Field* from);

FieldShape* getShapeByName(const char* name);

/** \brief count the number of nodes affecting an element type
  \param type select from apf::Mesh::Type */
int countElementNodes(FieldShape* s, int type);

/** \brief gets the xi coordinates for all the nodes
  \details order follows canonical notation. See tables
  apf::Mesh::tri_edge_verts, apf::Mesh::tet_edge_verts, and
  apf::Mesh::tet_tri_verts
  \param type select from apf::Mesh::Type */
void getElementNodeXis(FieldShape* s, int type,
    apf::NewArray<apf::Vector3>& xis);

/** \brief gets the xi coordinates for all the nodes
  \details order follows downward adjacency and global
  directions for the bounding entities. xi coordinates
  will be with respect to the entity e
  \param type select from apf::Mesh::Type */
void getElementNodeXis(FieldShape* s, Mesh* m, MeshEntity* e,
    apf::NewArray<apf::Vector3>& xis);

/** \brief Reparameterize from boundary entity to element
  \details This function converts a point in the local
  parametric space of a boundary mesh entity into the
  equivalent point in the local parametric space of
  an adjacent element.
  */
Vector3 boundaryToElementXi(
    Mesh* m,
    MeshEntity* boundary,
    MeshEntity* element,
    Vector3 const& xi);

}

#endif

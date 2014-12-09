/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFNUMBERING_H
#define APFNUMBERING_H

/** \file apfNumbering.h
  \brief Local and global numbering interface */

#include "apf.h"
#include "apfDynamicArray.h"
#include "apfMesh.h"

namespace apf {

/** \brief Global numberings use 64-bit integers */
typedef NumberingOf<long> GlobalNumbering;

/** \brief Create a Numbering of degrees of freedom of a Field.
  */
Numbering* createNumbering(Field* f);

/** \brief Create a generally-defined Numbering 
  * \details This numbering will be available via mesh->findNumbering, etc.
  *   The shape determines where the nodes are, and the component count
  *   determines how many integers there are per node.
  */
Numbering* createNumbering(
    Mesh* mesh,
    const char* name,
    FieldShape* shape,
    int components);

/** \brief Destroy a Numbering.
  */
void destroyNumbering(Numbering* n);

/** \brief Set the fixed/free status of a degree of freedom,
  * \details must be called prior to making any isFixed()
  *   calls on the same node/component.
  *
  * \param n the numbering object
  * \param e the mesh entity with which the node is associated
  * \param node the node number withing the mesh entity
  * \param component the component number within the nodal tensor
  */
void fix(Numbering* n, MeshEntity* e, int node, int component, bool fixed);

/** \brief Check whether a degree of freedom is fixed
  */
bool isFixed(Numbering* n, MeshEntity* e, int node, int component);

/** \brief Check whether a degree of freedom is numbered
  */
bool isNumbered(Numbering* n, MeshEntity* e, int node, int component);

/** \brief number a degree of freedom
  */
void number(Numbering* n, MeshEntity* e, int node, int component, int number);

/** \brief get a degree of freedom number
  */
int getNumber(Numbering* n, MeshEntity* e, int node, int component);

/** \brief get the field being numbered
  */
Field* getField(Numbering* n);

/** \brief get the FieldShape used by a Numbering */
FieldShape* getShape(Numbering* n);
/** \brief get the name of a Numbering */
const char* getName(Numbering* n);
/** \brief get the mesh associated with a Numbering */
Mesh* getMesh(Numbering* n);

/** \brief returns the node numbers of an element
  \details numbers are returned in the standard
           element node ordering for its shape functions */
void getElementNumbers(Numbering* n, MeshEntity* e, NewArray<int>& numbers);

/** \brief return the number of fixed degrees of freedom */ 
int countFixed(Numbering* n);

/** \brief numbers non-owned nodes with the values from their owners
   \details Works even if the non-owned nodes have no number currently
   assigned, after this they are all numbered
   \param shr if non-zero, use this Sharing model to determine ownership
              and copies, otherwise call apf::getSharing */
void synchronize(Numbering * n, Sharing* shr = 0);

/** \brief number the local owned entities of a given dimension */
Numbering* numberOwnedDimension(Mesh* mesh, const char* name, int dim);
/** \brief number the local elements */
Numbering* numberElements(Mesh* mesh, const char* name);
/** \brief number all local nodes
  \param s if non-zero, use nodes from this FieldShape, otherwise
           use the mesh's coordinate nodes */
Numbering* numberOverlapNodes(
		Mesh* mesh,
		const char* name,
		FieldShape* s = 0);
/** \brief number the local owned nodes
  \param s if non-zero, use nodes from this FieldShape, otherwise
           use the mesh's coordinate nodes
  \param shr if non-zero, use this Sharing to determine ownership,
             otherwise call apf::getSharing */
Numbering* numberOwnedNodes(
    Mesh* mesh,
    const char* name,
    FieldShape* s = 0,
    Sharing* shr = 0);
/** \brief count the number of nodes that have been numbered */
int countNodes(Numbering* n);
int countNodes(GlobalNumbering* n);

/** \brief Node identifier */
struct Node
{
  Node() {}
  Node(MeshEntity* e, int n):entity(e),node(n) {}
  /** \brief unique entity to which the node is associated */
  MeshEntity* entity;
  /** \brief which node on the entity it is */
  int node;
};

/** \brief get an array of numbered nodes
  \details the array is sorted by
  dimension first and then by mesh iterator order */
void getNodes(Numbering* n, DynamicArray<Node>& nodes);

/** \brief get nodes on the closure of a model entity
  \details
   all local nodes associated with mesh entities classified
   on the closure of the model entity are returned.
   this is useful for applying boundary conditions, and correctly
   handles cases when a part has, for example, vertices classified
   on the edge of a model face, but none of the triangles classified
   on the model face itself. */
void getNodesOnClosure(
    Mesh* m,
    ModelEntity* me,
    DynamicArray<Node>& on);

/** \brief create global numbering
   \details see apf::createNumbering.
   so far global numberings have one component */
GlobalNumbering* createGlobalNumbering(
    Mesh* mesh,
    const char* name,
    FieldShape* shape);
FieldShape* getShape(GlobalNumbering* n);
const char* getName(GlobalNumbering* n);
/** \brief get the mesh associated with a global numbering */
Mesh* getMesh(GlobalNumbering* n);
/** \brief assign a global number */
void number(GlobalNumbering* n, Node node, long number);
/** \brief get a global number */
long getNumber(GlobalNumbering* n, Node node);
/** \brief get an element's global node numbers */
int getElementNumbers(GlobalNumbering* n, MeshEntity* e,
    NewArray<long>& numbers);

/** \brief converts a local numbering into a global numbering.
  \details
   the original local numbering is destroyed.
   All local numbers are increased by an offset;
   the offset on part P is the sum of the numbered nodes
   on parts [0,P-1].

   this is done in O(log N) time in parallel, where N is the part count.
 
   this function has no intrinsic knowledge of ownership,
   it operates simply on nodes which have been explicitly numbered.
   the input to this function is usually a numbering produced by
   a numberOwned* function, and the result is a global numbering of
   all the owned nodes. subsequently calling synchronize on the global
   numbering completes the typical process which leaves all nodes
   (owned and not) with a global number attached.
 */
GlobalNumbering* makeGlobal(Numbering* n);

/** \brief see the Numbering equivalent and apf::makeGlobal */
void synchronize(GlobalNumbering* n, Sharing* shr = 0);

/** \brief destroy a global numbering */
void destroyGlobalNumbering(GlobalNumbering* n);

/** \brief see the Numbering equivalent */
void getNodes(GlobalNumbering* n, DynamicArray<Node>& nodes);

/** \brief Number by adjacency graph traversal
  \details a plain single-integer tag is used to
  number the vertices and elements of a mesh */
MeshTag* reorder(Mesh* mesh, const char* name);

void globalize(Numbering* n);

/** \brief number all components by simple iteration
 \todo name should be lower-case */
int NaiveOrder(Numbering * num);

/** \brief like apf::reorder, but numbers all free nodal components
 \todo name should be lower-case */
int AdjReorder(Numbering * num);

/** \brief add an offset to all free nodal component numbers
 \todo name should be lower-case */
void SetNumberingOffset(Numbering * num, int off);

}

#endif

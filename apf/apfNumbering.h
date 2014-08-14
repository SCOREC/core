/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFNUMBERING_H
#define APFNUMBERING_H

#include "apf.h"
#include "apfDynamicArray.h"
#include "apfMesh.h"

namespace apf {

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

FieldShape* getShape(Numbering* n);
const char* getName(Numbering* n);
Mesh* getMesh(Numbering* n);

/* returns the numbers of the nodes of an element, in the standard
   element node ordering for that element */
void getElementNumbers(Numbering* n, MeshEntity* e, NewArray<int>& numbers);
int countFixed(Numbering* n);

/* numbers non-owned nodes with the values from their owned copies on
   other parts. Works even if the non-owned nodes have no number currently
   assigned, after this they are all numbered */
void synchronize(Numbering * n, Sharing* shr = 0);

/* creates a local numbering of the owned entities of a given dimension */
Numbering* numberOwnedDimension(Mesh* mesh, const char* name, int dim);
/* creates a local element numbering */
Numbering* numberElements(Mesh* mesh, const char* name);
/* creates a local numbering of all nodes (owned or not)
   in the mesh's default field shape.
   usually this is the vertices. */
Numbering* numberOverlapNodes(
		Mesh* mesh,
		const char* name,
		FieldShape* s = 0);
/* creates a local numbering of the owned nodes in the
   mesh's default field shape. */
Numbering* numberOwnedNodes(
    Mesh* mesh,
    const char* name,
    FieldShape* s = 0,
    Sharing* shr = 0);
/* count the number of nodes that have been explicitly numbered */
int countNodes(Numbering* n);

struct Node
{
  Node() {}
  Node(MeshEntity* e, int n):entity(e),node(n) {}
  MeshEntity* entity;
  int node;
};

/* returns an array of all local nodes (owned or not) with
   an explicit number assignment. the array is sorted by
   dimension first and then by mesh iterator order */
void getNodes(Numbering* n, DynamicArray<Node>& nodes);

/* returns an array of all nodes associated with entities which
   are classified on the closure of a model face.
   this is useful for applying boundary conditions, and correctly
   handles cases when a part has, for example, vertices classified
   on the edge of a model face, but none of the triangles classified
   on the model face itself. */
void getNodesOnClosure(
    Mesh* m,
    ModelEntity* me,
    DynamicArray<Node>& on);

/* creates an empty global numbering, see createNumbering.
   so far global numberings have one component */
GlobalNumbering* createGlobalNumbering(
    Mesh* mesh,
    const char* name,
    FieldShape* shape);
/* see the Numbering equivalents of these functions */
Mesh* getMesh(GlobalNumbering* n);
void number(GlobalNumbering* n, Node node, long number);
long getNumber(GlobalNumbering* n, Node node);
int getElementNumbers(GlobalNumbering* n, MeshEntity* e,
    NewArray<long>& numbers);

/* converts a local numbering into a global numbering.
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

/* see synchronize(Numbering*) and makeGlobal */
void synchronize(GlobalNumbering* n, Sharing* shr = 0);

void destroyGlobalNumbering(GlobalNumbering* n);

/* see the Numbering equivalent */
void getNodes(GlobalNumbering* n, DynamicArray<Node>& nodes);

}

#endif

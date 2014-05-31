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

Numbering* createNumbering(
    Mesh* mesh,
    const char* name,
    FieldShape* shape,
    int components);

/** \brief Destroy a Numbering.
  */
void destroyNumbering(Numbering* n);

/** \brief Set the fixed/free status of a degree of freedom, must be called prior to making any isFixed() calls on the same node/component.
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

void getElementNumbers(Numbering* n, MeshEntity* e, NewArray<int>& numbers);
int countFixed(Numbering* n);

void synchronize(Numbering * n);

Numbering* numberElements(Mesh* mesh, const char* name);
Numbering* numberOverlapNodes(
		Mesh* mesh,
		const char* name,
		FieldShape* s = 0);
Numbering* numberOwnedNodes(
    Mesh* mesh,
    const char* name,
    FieldShape* s = 0);
int countNodes(Numbering* n);

struct Node
{
  Node() {}
  Node(MeshEntity* e, int n):entity(e),node(n) {}
  MeshEntity* entity;
  int node;
};

void getNodes(Numbering* n, DynamicArray<Node>& nodes);

void getNodesOnClosure(
    Mesh* m,
    ModelEntity* me,
    DynamicArray<Node>& on);

GlobalNumbering* createGlobalNumbering(
    Mesh* mesh,
    const char* name,
    FieldShape* shape);
Mesh* getMesh(GlobalNumbering* n);
void number(GlobalNumbering* n, Node node, long number);
long getNumber(GlobalNumbering* n, Node node);
int getElementNumbers(GlobalNumbering* n, MeshEntity* e, NewArray<long>& numbers);
GlobalNumbering* makeGlobal(Numbering* n);
void synchronize(GlobalNumbering* n);
void destroyGlobalNumbering(GlobalNumbering* n);
void getNodes(GlobalNumbering* n, DynamicArray<Node>& nodes);

/* a bit of backward compatibility, main user is Epetra Albany */
void globalize(Numbering* n);

}

#endif

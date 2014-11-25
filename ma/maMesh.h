/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_MESH_H
#define MA_MESH_H

/** \file maMesh.h
  \brief mesh functions for MeshAdapt */

#include <apfMesh2.h>
#include <apfMatrix.h>
#include <set>

namespace ma {

/** \brief convenient vector name */
typedef apf::Vector3 Vector;
/** \brief convenient matrix name */
typedef apf::Matrix3x3 Matrix;
/** \brief convenient mesh name */
typedef apf::Mesh2 Mesh;
/** \brief convenient mesh entity name */
typedef apf::MeshEntity Entity;
/** \brief convenient mesh iterator name */
typedef apf::MeshIterator Iterator;
/** \brief convenient mesh tag name */
typedef apf::MeshTag Tag;
/** \brief convenient mesh entity array name */
typedef apf::DynamicArray<Entity*> EntityArray;
/** \brief convenient mesh entity set name */
typedef std::set<Entity*> EntitySet;
/** \brief convenient mesh entity upward adjacencies name */
typedef EntityArray Upward;
/** \brief convenient mesh entity downward adjacencies name */
typedef apf::Downward Downward;
/** \brief convenient geometric model entity name */
typedef apf::ModelEntity Model;

enum MeshEntityType
{
  VERT = apf::Mesh::VERTEX,      //0
  EDGE = apf::Mesh::EDGE,        //1
  TRI = apf::Mesh::TRIANGLE,     //2
  QUAD = apf::Mesh::QUAD,        //3
  TET = apf::Mesh::TET,          //4
  HEX = apf::Mesh::HEX,          //5
  PRISM = apf::Mesh::PRISM,      //6
  PYRAMID = apf::Mesh::PYRAMID,  //7
  TYPES = apf::Mesh::TYPES       //8
};

/** \brief get vertex spatial coordinates */
Vector getPosition(Mesh* m, Entity* vertex);

/** \brief convenient remote copies name */
typedef apf::Copies Remotes;
/** \brief part id set name */
typedef apf::Parts Parts;

void rotateTri(Entity** iv, int n, Entity** ov);
void rotateQuad(Entity** iv, int n, Entity** ov);
void rotateTet(Entity** iv, int n, Entity** ov);
void rotatePrism(Entity** iv, int n, Entity** ov);
void rotatePyramid(Entity** iv, int n, Entity** ov);
void rotateEntity(int type, Entity** iv, int n, Entity** ov);
void rotateEntity(apf::Mesh* m, Entity* e, int n, Entity** v);

int findTetRotation(Mesh* m, Entity* tet, Entity** v);
void unrotateTetXi(Vector& xi, int rotation);

void rotateOct(Entity** iv, int n, Entity** ov);

int getDownIndex(Mesh* m, Entity* e, Entity* de);
Entity* getTriEdgeOppositeVert(Mesh* m, Entity* tri, Entity* v);
Entity* getTriVertOppositeEdge(Mesh* m, Entity* tri, Entity* v);
Entity* getTetVertOppositeTri(Mesh* m, Entity* tet, Entity* tri);
Entity* getQuadEdgeOppositeEdge(Mesh* m, Entity* q, Entity* e);

Entity* findTetByTwoTris(Mesh* m, Entity** tris);

/** \brief rebuild an element with one vertex being different
  \details uses the original
  to reconstruct geometric classification */
Entity* rebuildElement(
    Mesh* m,
    Entity* original,
    Entity* oldVert,
    Entity* newVert,
    apf::BuildCallback* cb);

bool isInClosure(Mesh* m, Entity* parent, Entity* e);

void getBoundingBox(Mesh* m, Vector& lower, Vector& upper);
Vector getCentroid(Mesh* m);

void ensureParallelConsistency(Mesh* m);

Entity* findTriFromVerts(Mesh* m, Entity** v);

double measure(Mesh* m, Entity* e);

bool isOnModelEdge(Mesh* m, Entity* e);
bool isOnModelFace(Mesh* m, Entity* e);

Vector getTriNormal(Mesh* m, Entity** v);
Vector getTriNormal(Mesh* m, Entity* e);
bool isTwoTriAngleAcute(Mesh* m, Entity** va, Entity** vb);
bool isTwoTriAngleAcute(Mesh* m, Entity* a, Entity* b);

/** \brief Computes the insphere radius of an element
  * \todo currently only implemented for tets
  */
double getInsphere(Mesh* m, Entity* e);

double getAverageElementSize(Mesh* m);
double getMinimumElementSize(Mesh* m);

void getFaceEdgesAndDirections(
    Mesh* m,
    Entity* face,
    Entity** edges,
    int* directions);

int getFaceEdgeDirection(
    Mesh* m,
    Entity* face,
    Entity* edge);

Entity* findEdge(Mesh* m, Entity* v0, Entity* v1);
bool edgeExists(Mesh* m, Entity* v0, Entity* v1);

/* returns true iff the direction of the edge is along
   the direction of the triangle's vertex (and edge) ordering */
bool isTriEdgeAligned(Mesh* m, Entity* tri, Entity* edge);

}

#endif

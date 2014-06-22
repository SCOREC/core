/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_MESH_H
#define MA_MESH_H

#include <apfMesh2.h>
#include <apfMatrix.h>
#include <set>

namespace ma {

typedef apf::Vector3 Vector;
typedef apf::Matrix3x3 Matrix;
typedef apf::Mesh2 Mesh;
typedef apf::MeshEntity Entity;
typedef apf::MeshIterator Iterator;
typedef apf::MeshTag Tag;
typedef apf::DynamicArray<Entity*> EntityArray;
typedef std::set<Entity*> EntitySet;
typedef EntityArray Upward;
typedef apf::Downward Downward;
typedef apf::ModelEntity Model;

enum {
  VERT = apf::Mesh::VERTEX,
  EDGE = apf::Mesh::EDGE,
  TRI = apf::Mesh::TRIANGLE,
  QUAD = apf::Mesh::QUAD,
  TET = apf::Mesh::TET,
  HEX = apf::Mesh::HEX,
  PRISM = apf::Mesh::PRISM,
  PYRAMID = apf::Mesh::PYRAMID,
  TYPES = apf::Mesh::TYPES
};

Vector getPosition(Mesh* m, Entity* vertex);

typedef apf::Copies Remotes;
typedef apf::Parts Parts;

void addRemote(Mesh* m, Entity* e, int p, Entity* r);

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

/* averages the real coordinates of n vertices together */
Vector averagePositions(Mesh* m, Entity** v, int n);

int getDownIndex(Mesh* m, Entity* e, Entity* de);
Entity* getTriEdgeOppositeVert(Mesh* m, Entity* tri, Entity* v);
Entity* getTriVertOppositeEdge(Mesh* m, Entity* tri, Entity* v);
Entity* getTetVertOppositeTri(Mesh* m, Entity* tet, Entity* tri);
Entity* getQuadEdgeOppositeEdge(Mesh* m, Entity* q, Entity* e);

Entity* findTetByTwoTris(Mesh* m, Entity** tris);

/* constructs a new element using the original
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

Vector getTriNormal(Mesh* m, Vector* x);
bool isTwoTriAngleAcute(Mesh* m, Entity** va, Entity** vb);
bool isTwoTriAngleAcute(Mesh* m, Entity* a, Entity* b);

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

}

#endif

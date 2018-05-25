/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_SNAP_H
#define MA_SNAP_H

#include "maMesh.h"

namespace ma {

class Adapt;

void snap(Adapt* a);
void visualizeGeometricInfo(Mesh* m, const char* name);

long snapTaggedVerts(Adapt* a, Tag* snapTag);

void interpolateParametricCoordinates(
    apf::Mesh* m,
    Model* g,
    double t,
    Vector const& a,
    Vector const& b,
    Vector& p);
void transferParametricOnEdgeSplit(
    Mesh* m,
    Entity* e,
    double t,
    Vector& p);
void transferParametricOnTriSplit(
    Mesh* m,
    Entity* face,
    const Vector& xi,
    Vector& param);
void transferParametricOnQuadSplit(
    Mesh* m,
    Entity* quad,
    Entity* v01,
    Entity* v32,
    double y,
    Vector& p);

void getClosestPointParametricCoordinates(
    apf::Mesh* m,
    Model* g,
    double t,
    Vector const& a,
    Vector const& b,
    Vector& p);
void transferToClosestPointOnEdgeSplit(
    Mesh* m,
    Entity* e,
    double t,
    Vector& p);
void transferToClosestPointOnTriSplit(
    Mesh* m,
    Entity* face,
    const Vector& xi,
    Vector& param);
void transferToClosestPointOnQuadSplit(
    Mesh* m,
    Entity* quad,
    Entity* v01,
    Entity* v32,
    double y,
    Vector& p);

}

#endif

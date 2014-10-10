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

void transferParametricOnEdgeSplit(
    Mesh* m,
    Entity* e,
    double t,
    Vector& p);
void transferParametricOnQuadSplit(
    Mesh* m,
    Entity* quad,
    Entity* v01,
    Entity* v32,
    double y,
    Vector& p);

}

#endif

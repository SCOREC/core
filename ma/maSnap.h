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

class Refine;

void transferParametricOnEdgeSplit(
    Mesh* m,
    Entity* e,
    double t,
    Vector& p);
void snap(Refine* r);
void visualizeGeometricInfo(Mesh* m, const char* name);

}

#endif

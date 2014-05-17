/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/

#ifndef MA_MAP_H
#define MA_MAP_H

#include "maAffine.h"
#include "maMesh.h"

namespace ma {

void getVertPoints(apf::Mesh* m, Entity* e, Vector* p);

Affine getMap(apf::Mesh* m, Entity* e);

double getInsideness(apf::Mesh* m, Entity* e, Vector const& xi);

}

#endif


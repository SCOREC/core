/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_EDGESWAP_H
#define MA_EDGESWAP_H

#include "maMesh.h"

namespace ma {

class Adapt;

class EdgeSwap
{
  public:
    virtual ~EdgeSwap() {};
    virtual bool run(Entity* edge) = 0;
};

EdgeSwap* makeEdgeSwap(Adapt* a);

}

#endif

/****************************************************************************** 

  Copyright (c) 2004-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_MATCH_H
#define MA_MATCH_H

#include "maMesh.h"

namespace ma {

class Refine;

void matchNewVerts(Refine* r);
void matchNewElements(Refine* r);

}

#endif

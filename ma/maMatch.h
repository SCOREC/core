/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_MATCH_H
#define MA_MATCH_H

namespace ma {

class Refine;
class Adapt;

void matchNewVerts(Refine* r);
void matchNewElements(Refine* r);

void preventMatchedCavityMods(Adapt* a);

}

#endif

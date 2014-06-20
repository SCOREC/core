/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_SNAPPER_H
#define MA_SNAPPER_H

#include "maCollapse.h"

namespace apf {
class CavityOp;
}

namespace ma {

/* tries to make room for a vertex to snap into
   a mesh by collapsing edges in the direction
   of desired snapping */

class Snapper
{
  public:
    Snapper(Adapt* a, Tag* st, bool is);
    bool setVert(Entity* v, apf::CavityOp* o);
    bool run();
    bool dug;
  private:
    Adapt* adapter;
    Tag* snapTag;
    Entity* vert;
    Collapse collapse;
    bool isSimple;
};

}

#endif

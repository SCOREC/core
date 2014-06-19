/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_DIGGER_H
#define MA_DIGGER_H

#include "maCollapse.h"

namespace apf {
class CavityOp;
}

namespace ma {

/* tries to make room for a vertex to snap into
   a mesh by collapsing edges in the direction
   of desired snapping */

class Digger
{
  public:
    Digger(Adapt* a, Tag* st);
    /* this will request the vertex and all
       vertices adjacent to it through an edge,
       i.e. two layers of elements. That allows
       us to try all possible collapses in that cavity.*/
    bool setVert(Entity* v, apf::CavityOp* o);
    bool run();
  private:
    bool tryToCollapse(Entity* e);
    Adapt* adapter;
    Mesh* mesh;
    Tag* snapTag;
    Entity* vert;
    EntityArray edges;
    Collapse collapse;
};

}

#endif

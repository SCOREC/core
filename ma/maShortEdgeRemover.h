/****************************************************************************** 

  Copyright (c) 2004-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_SHORTEDGEREMOVER_H
#define MA_SHORTEDGEREMOVER_H

/* short edge removal is related
   to edge collapsing but not identical.
   Here we try to eliminate the short edge
   by any means necessary, typically because
   collapsing it has failed.
   
   We begin with an approach suggested by
   Jie Wan: try to collapse adjacent edges. */

/* note - this algorithm is quite agressive, and requires two
   layers of elements around the target vertices.
   It should be run sparingly. */

#include "maVertRemover.h"

namespace ma {

class ShortEdgeRemover
{
  public:
    ShortEdgeRemover(Adapt* a);
    void setEdge(Entity* e);
    bool requestLocality(apf::CavityOp* o);
    void findEdges();
    void getOldElements(EntityArray& oldElements);
    bool didImproveQuality();
    bool tryToCollapse(Entity* e, Entity* v);
    bool tryToRemoveVert(int vi);
    bool run();
  private:
    Adapt* adapter;
    Mesh* mesh;
    Entity* edge;
    VertRemover vertRemovers[2];
};

}

#endif

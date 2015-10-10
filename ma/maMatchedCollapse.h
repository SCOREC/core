/****************************************************************************** 

  Copyright 2015 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_MATCHED_COLLAPSE_H
#define MA_MATCHED_COLLAPSE_H

#include "maCollapse.h"

namespace ma {

struct MatchedCollapse
{
  MatchedCollapse(Adapt* a);
  void setEdge(Entity* e);
  bool requestLocality(apf::CavityOp* o);
  void setEdges();
  bool checkTopo();
  Adapt* adapt;
  Mesh* mesh;
  apf::Sharing* sharing;
  apf::DynamicArray<Collapse> collapses;
};

}

#endif

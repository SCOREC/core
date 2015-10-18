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

struct Rebuild {
  Rebuild(Entity* a, Entity* b);
  Entity* e;
  Entity* original;
  bool operator<(Rebuild const&other) const;
  bool operator==(Rebuild const&other) const;
};

struct Rebuilds : public RebuildCallback {
  Rebuilds(Mesh* m);
  virtual void rebuilt(Entity* e, Entity* original);
  void reset();
  void match(apf::Sharing* sh,
      apf::DynamicArray<ma::Collapse>& collapses);
  Mesh* mesh;
  std::vector<Rebuild> v;
};

struct MatchedCollapse
{
  MatchedCollapse(Adapt* a);
  void setEdge(Entity* e);
  bool requestLocality(apf::CavityOp* o);
  void setEdges();
  bool checkTopo();
  void unmark();
  void cancel();
  bool tryThisDirection(double qualityToBeat);
  bool tryBothDirections(double qualityToBeat);
  void destroyOldElements();
  Adapt* adapt;
  Mesh* mesh;
  apf::Sharing* sharing;
  apf::DynamicArray<Collapse> collapses;
  Rebuilds rebuilds;
};

}

#endif

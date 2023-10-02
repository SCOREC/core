/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_MATCHED_SNAPPER_H
#define MA_MATCHED_SNAPPER_H

#include "maSnapper.h"

namespace apf {
class CavityOp;
}

namespace ma {

class MatchedSnapper
{
  public:
    MatchedSnapper(Adapt* a, Tag* st, bool is);
    ~MatchedSnapper();
    void setVert(Entity* v);
    bool requestLocality(apf::CavityOp* o);
    void setVerts();
    bool trySnaps();
    void cancelSnaps();
  private:
    Adapt* adapter;
    Tag* snapTag;
    Entity* vert;
    apf::DynamicArray<Snapper*> snappers;
    apf::DynamicArray<Vector> locations;
    bool isSimple;
    apf::Sharing* sharing;
};

}

#endif

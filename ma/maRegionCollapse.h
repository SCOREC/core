/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_REGION_COLLAPSE_H
#define MA_REGION_COLLAPSE_H

#include "maAdapt.h"

namespace apf {
class CavityOp;
}

namespace ma {

class Adapt;

class RegionCollapse
{
  public:
    void Init(Adapt* a, double fa);
    bool requestLocality(apf::CavityOp* o);
    bool setRegion(Entity* r);
    bool checkGeom();
    bool checkTopo();
    void apply();
    void unmark();
    Adapt* adapt;
    Entity* region;
    Entity* reclassifyEdge;
    int numBdryFaces;
    Entity* faces[4];
    double flatAngle;
};

bool setupRegionCollapse(RegionCollapse& rcollapse, Entity* region);

}

#endif

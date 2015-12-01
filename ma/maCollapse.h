/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_COLLAPSE_H
#define MA_COLLAPSE_H

#include "maAdapt.h"

namespace apf {
class CavityOp;
}

namespace ma {

class Adapt;

class Collapse
{
  public:
    void Init(Adapt* a);
    bool requestLocality(apf::CavityOp* o);
    void destroyOldElements();
    void destroyNewElements();
    bool setEdge(Entity* e);
    bool checkClass();
    bool checkTopo();
    void unmark();
    void setVerts();
    void computeElementSets();
    void rebuildElements();
    bool isGood2DMesh();
    void cancel();
    bool tryThisDirection(double qualityToBeat);
    bool tryThisDirectionNoCancel(double qualityToBeat);
    bool tryBothDirections(double qualityToBeat);
    void getOldElements(EntityArray& oldElements);
    double getOldQuality();
    Adapt* adapt;
    Entity* edge; 
    Entity* vertToCollapse;
    Entity* vertToKeep;
    EntitySet elementsToCollapse;
    EntitySet elementsToKeep;
    EntityArray newElements;
    Cavity cavity;
    RebuildCallback* rebuildCallback;
};

bool checkEdgeCollapseTopology(Adapt* a, Entity* edge);
bool isRequiredForAnEdgeCollapse(Adapt* adapt, Entity* vertex);
bool setupCollapse(Collapse& collapse, Entity* edge, Entity* vert);

}

#endif

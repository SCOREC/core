/****************************************************************************** 

  Copyright 2017 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_FACESPLITCOLLAPSE_H
#define MA_FACESPLITCOLLAPSE_H

#include "maFaceSplit.h"
#include "maCollapse.h"

namespace ma {

class FaceSplitCollapse
{
  class IgnoringCollapse : public Collapse
  {
  public:
    EntitySet elementsToIgnore;
    virtual void computeElementSets();
  };

  public:
    FaceSplitCollapse(Adapt* a);
    void getNewElements(EntityArray& e);
    bool tryBothCollapses(Entity* e);
    void accept();
    bool run(Entity* face, Entity* tet);
    Adapt* getAdapt();
  private:
    FaceSplit faceSplit;
    IgnoringCollapse collapse;
    double oldQuality;
};

}

#endif

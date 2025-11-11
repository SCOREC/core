/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_DOUBLESPLITCOLLAPSE_H
#define MA_DOUBLESPLITCOLLAPSE_H

#include "maSplits.h"
#include "maCollapse.h"

namespace ma {

class DoubleSplitCollapse
{
  public:
    DoubleSplitCollapse(Adapt* a);
    void getNewElements(EntityArray& e);
    bool didImproveQuality();
    bool tryThisCollapse();
    bool tryBothCollapses(Entity* e);
    void accept();
    /* Quality is optional since during fix shape we don't need a target quality as
    the quality improves. However we need a target quality during snapping since
    we don't care if the quality gets worse. */
    bool run(Entity** edges, double quality = -1);
    Adapt* getAdapt();
  private:
    Splits splits;
    Collapse collapse;
    double oldQuality;
};

}

#endif

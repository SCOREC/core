/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_SINGLESPLITCOLLAPSE_H
#define MA_SINGLESPLITCOLLAPSE_H

#include "maSplits.h"
#include "maCollapse.h"

namespace ma {

class SingleSplitCollapse
{
  class IgnoringCollapse : public Collapse
  {
  public:
    EntitySet elementsToIgnore;
    virtual void computeElementSets();
    virtual bool setEdge(Entity* e);
  };

public:
  SingleSplitCollapse(Adapt* a);
  void getNewElements(EntityArray& e);
  bool didImproveQuality();
  bool tryThisCollapse();
  bool tryBothCollapses(Entity* e);
  void accept();
  /* Quality is optional since during fix shape we don't need a target quality as
  the quality improves. However we need a target quality during snapping since
  we don't care if the quality gets worse. */
  bool run(Entity* edge, Entity* vert, double quality=-1);
  Adapt* getAdapt();
private:
  Entity *oldVert, *oldEdge;
  Splits splits;
  IgnoringCollapse collapse;
  double oldQuality;
};

}

#endif

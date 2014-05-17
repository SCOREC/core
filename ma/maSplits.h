/****************************************************************************** 

  Copyright (c) 2004-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_SPLIT_H
#define MA_SPLIT_H

#include "maRefine.h"

namespace ma {

class Splits
{
  public:
    Splits(Adapt* a);
    bool setEdges(Entity** e, int n);
    void makeNewElements();
    void cancel();
    void transfer();
    void destroyOldElements();
    Entity* getSplitVert(int i);
    EntityArray& getTets() {return refiner->toSplit[3];}
    Adapt* getAdapt() {return refiner->adapt;}
  private:
    Refine* refiner;
};

}

#endif

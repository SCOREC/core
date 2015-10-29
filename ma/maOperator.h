/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_MODIFIER_H
#define MA_MODIFIER_H

#include "maAdapt.h"
#include <apfCavityOp.h>

namespace ma {

class Operator
{
  public:
    virtual ~Operator();
    virtual int getTargetDimension() = 0;
    virtual bool shouldApply(Entity* e) = 0;
    virtual bool requestLocality(apf::CavityOp* o) = 0;
    virtual void apply() = 0;
};

void applyOperator(Adapt* a, Operator* o, bool matched = false);

}

#endif

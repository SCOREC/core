/****************************************************************************** 

  Copyright (c) 2004-2014 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maOperator.h"
#include "maAdapt.h"

namespace ma {

class CollectiveOperation : public apf::CavityOp, public DeleteCallback
{
  public:
    CollectiveOperation(Adapt* a, Operator* o):
      apf::CavityOp(a->mesh,true),
      DeleteCallback(a)
    {
      op = o;
    }
    Outcome setEntity(Entity* e)
    {
      if ( ! op->shouldApply(e))
        return SKIP;
      if ( ! op->requestLocality(this))
        return REQUEST;
      return OK;
    }
    void apply()
    {
      op->apply();
    }
    void call(Entity* e)
    {
      this->preDeletion(e);
    }
  private:
    Operator* op;
};

Operator::~Operator() {}

void applyOperator(Adapt* a, Operator* o)
{
  CollectiveOperation op(a,o);
  op.applyToDimension(o->getTargetDimension());
}

}

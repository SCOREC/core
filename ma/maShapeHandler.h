#ifndef MA_SHAPEHANDLER_H
#define MA_SHAPEHANDLER_H

#include "maSolutionTransfer.h"

namespace ma {

class Adapt;

class ShapeHandler : public SolutionTransfer
{
  public:
    virtual double getQuality(Entity* e) = 0;
};

ShapeHandler* getShapeHandler(Adapt* a);

}

#endif

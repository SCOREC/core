/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFVECTORFIELD_H
#define APFVECTORFIELD_H

#include "apfFieldOf.h"
#include "apf.h"

namespace apf {

class VectorField : public FieldOf<Vector3>
{
  public:
    virtual ~VectorField() {}
    virtual Element* getElement(VectorElement* e);
    virtual int getValueType() const {return VECTOR;}
    virtual int getShapeType() const {return SCALAR;}
    virtual int countComponents() const;
};

}//namespace apf

#endif

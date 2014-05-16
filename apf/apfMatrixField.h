/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFMATRIXFIELD_H
#define APFMATRIXFIELD_H

#include "apfFieldOf.h"
#include "apf.h"

namespace apf {

class MatrixField : public FieldOf<Matrix3x3>
{
  public:
    virtual ~MatrixField() {}
    virtual Element* getElement(VectorElement* e);
    virtual int getValueType() const {return MATRIX;}
    virtual int countComponents() const;
};

}//namespace apf

#endif

/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFMATRIXELEMENT_H
#define APFMATRIXELEMENT_H

#include "apfElementOf.h"

namespace apf {

class MatrixField;

class MatrixElement : public ElementOf<Matrix3x3>
{
  public:
    MatrixElement(MatrixField* f, VectorElement* e);
    virtual ~MatrixElement();
};

}//namespace apf

#endif

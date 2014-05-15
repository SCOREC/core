/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFSCALARFIELD_H
#define APFSCALARFIELD_H

#include "apfFieldOf.h"
#include "apf.h"

namespace apf {

class ScalarField : public FieldOf<double>
{
  public:
    virtual ~ScalarField() {}
    virtual Element* getElement(VectorElement* e);
    virtual int getValueType() const {return SCALAR;}
    virtual int countComponents() const;
};

}//namespace apf

#endif

/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFPACKEDFIELD_H
#define APFPACKEDFIELD_H

#include "apfField.h"
#include "apf.h"

namespace apf {

class PackedField : public Field
{
  public:
    PackedField(int c):components(c) {}
    virtual ~PackedField() {}
    virtual Element* getElement(VectorElement* e);
    virtual int getValueType() const {return PACKED;}
    virtual int countComponents() const {return components;}
    virtual void project(Field* from);
    virtual void axpy(double a, Field* x);
  private:
    int components;
};

}//namespace apf

#endif


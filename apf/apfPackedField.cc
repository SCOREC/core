/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfPackedField.h"
#include "apfElement.h"

namespace apf {

Element* PackedField::getElement(VectorElement* e)
{
  return new Element(this,e);
}

void PackedField::project(Field* from)
{
  fail("PackedField::project unimplemented");
}

void PackedField::axpy(double a, Field* x)
{
  fail("PackedField::axpy unimplemented");
}

}//namespace apf


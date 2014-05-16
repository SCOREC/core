/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfVectorField.h"
#include "apfVectorElement.h"

namespace apf {

Element* VectorField::getElement(VectorElement* e)
{
  return new VectorElement(this,e);
}

int VectorField::countComponents() const
{
  return 3;
}

}//namespace apf

/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfMixedVectorField.h"
#include "apfMixedVectorElement.h"
#include "apfVectorElement.h"

namespace apf {

Element* MixedVectorField::getElement(VectorElement* e)
{
  return new MixedVectorElement(this,e);
}

int MixedVectorField::countComponents() const
{
  return 1;
}

}//namespace apf

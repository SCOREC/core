/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfScalarField.h"
#include "apfScalarElement.h"

namespace apf {

Element* ScalarField::getElement(VectorElement* e)
{
  return new ScalarElement(this,e);
}

int ScalarField::countComponents() const
{
  return 1;
}

}//namespace apf

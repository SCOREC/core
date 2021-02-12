/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfMatrixField.h"
#include "apfMatrixElement.h"

namespace apf {

Element* MatrixField::getElement(VectorElement* e)
{
  return new MatrixElement(this,e);
}

int MatrixField::countComponents() const
{
  return 9;
}

}//namespace apf

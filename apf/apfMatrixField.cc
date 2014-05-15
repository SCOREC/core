/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfMatrixField.h"
#include "apfMatrixElement.h"

namespace apf {

MatrixElement::MatrixElement(MatrixField* f, MeshElement* e):
    ElementOf<Matrix3x3>(f,e)
{
}

MatrixElement::~MatrixElement()
{
}

Element* MatrixField::getElement(VectorElement* e)
{
  return new MatrixElement(this,e);
}

int MatrixField::countComponents() const
{
  return 9;
}

}//namespace apf

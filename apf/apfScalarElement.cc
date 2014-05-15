/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfScalarElement.h"
#include "apfScalarField.h"

namespace apf {

ScalarElement::ScalarElement(ScalarField* f, VectorElement* e):
  ElementOf<double>(f,e)
{
}

void ScalarElement::grad(Vector3 const& local, Vector3& g)
{
  NewArray<Vector3> globalGradients;
  getGlobalGradients(local,globalGradients);
  double* nodeValues = getNodeValues();
  g = globalGradients[0] * nodeValues[0];
  for (int i=1; i < nen; ++i)
    g = g + globalGradients[i] * nodeValues[i];
}

}//namespace apf

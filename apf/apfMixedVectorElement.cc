/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <pcu_util.h>

#include "apfVectorElement.h"
#include "apfMixedVectorElement.h"
#include "apfMixedVectorField.h"

namespace apf {

MixedVectorElement::MixedVectorElement(MixedVectorField* f, MeshEntity* e):
  ElementOf<Vector3, double>(f,e)
{
}

MixedVectorElement::MixedVectorElement(MixedVectorField* f, VectorElement* p):
  ElementOf<Vector3, double>(f,p)
{
}

void MixedVectorElement::curl(Vector3 const& xi, Vector3& c)
{
  PCU_ALWAYS_ASSERT_VERBOSE( field->getShape()->isVectorShape(),
      "Not applicable for non-vector shape functions!");
    NewArray<Vector3> curlShapeValues;
    curlShapeValues.allocate(nen);
    getCurlShapeValues(this, xi, curlShapeValues);
    for (int ci = 0; ci < 3; ci++)
      c[ci] = 0.;
    for (int ni = 0; ni < nen; ni++)
      for (int ci = 0; ci < 3; ci++)
      	c[ci] += nodeData[ni] * curlShapeValues[ni][ci];
}

}//namespace apf

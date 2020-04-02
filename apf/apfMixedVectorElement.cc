/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

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
  //TODO to be completed using the curl of nedelec shapes and the nodeData
  (void) xi;
  (void) c;
}

}//namespace apf

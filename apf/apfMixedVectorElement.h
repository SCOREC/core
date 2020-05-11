/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFMIXEDVECTORELEMENT_H
#define APFMIXEDVECTORELEMENT_H

#include "apfElementOf.h"
#include "apfMatrix.h"

namespace apf {

class MixedVectorField;

/* Fields with vector shapes are a bit peculiar, in that
 * the shapes functions are vectors but the dof holders are
 * scalars. Hence the need for this Mixed class. An example of
 * such fields are Nedelec fields.
 */
class MixedVectorElement : public ElementOf<Vector3, double>
{
  public:
    MixedVectorElement(MixedVectorField* f, MeshEntity* e);
    MixedVectorElement(MixedVectorField* f, VectorElement* p);
    virtual ~MixedVectorElement() {}
    void curl(Vector3 const& xi, Vector3& c);
};

}//namespace apf

#endif

/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFSCALARELEMENT_H
#define APFSCALARELEMENT_H

#include "apfElementOf.h"

namespace apf {

class ScalarField;

class ScalarElement : public ElementOf<double>
{
  public:
    ScalarElement(ScalarField* f, VectorElement* e);
    virtual ~ScalarElement() {}
    void grad(Vector3 const& xi, Vector3& g);
};

}//namespace apf

#endif

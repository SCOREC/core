/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFMIXEDVECTORFIELD_H
#define APFMIXEDVECTORFIELD_H

#include "apfFieldOf.h"
#include "apf.h"

namespace apf {

/* This is used for fields with vector shape functions (e.g. Nedelec)
 * They are special in the sense that the dofs (or nodal values)
 * are scalar but the shape functions are vectors!
 */
class MixedVectorField : public FieldOf<double>
{
  public:
    virtual ~MixedVectorField() {}
    virtual Element* getElement(VectorElement* e);
    virtual int getValueType() const {return SCALAR;}
    virtual int getShapeType() const {return VECTOR;}
    virtual int countComponents() const;
};

}//namespace apf

#endif

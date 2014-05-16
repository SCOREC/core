/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFVECTORELEMENT_H
#define APFVECTORELEMENT_H

#include "apfElementOf.h"
#include "apfMatrix.h"

namespace apf {

class VectorField;

class VectorElement : public ElementOf<Vector3>
{
  public:
    VectorElement(VectorField* f, MeshEntity* e);
    VectorElement(VectorField* f, VectorElement* p);
    virtual ~VectorElement() {}
    double div(Vector3 const& xi);
    void grad(Vector3 const& xi, Matrix3x3& g);
    void curl(Vector3 const& xi, Vector3& c);
    void getJacobian(Vector3 const& xi, Matrix3x3& J);
    double getDV(Vector3 const& xi);
    void gradHelper(NewArray<Vector3>& nodalGradients, Matrix3x3& g);
};

double getJacobianDeterminant(Matrix3x3 const& J, int dimension);

}//namespace apf

#endif

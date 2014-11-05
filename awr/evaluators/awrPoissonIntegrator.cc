/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrPoissonIntegrator.h"

namespace awr {

PoissonIntegrator::PoissonIntegrator(apf::Field* f, int o) :
  apf::Integrator(o),
  f_(f) 
{
  shape_ = apf::getShape(f_);
  numComponents_ = apf::countComponents(f_);
  assert(numComponents_ == 1);
}

void PoissonIntegrator::inElement(apf::MeshElement* me)
{
  me_ = me;
  e_ = apf::createElement(f_,me_);
  numDims_ = apf::getDimension(me_);
  numNodes_ = apf::countNodes(e_);
  numLocalEqs_ = numNodes_*numComponents_;
  Ke_.setSize(numLocalEqs_,numLocalEqs_);
  for (int i=0; i < numLocalEqs_; ++i)
  for (int j=0; j < numLocalEqs_; ++j)
    Ke_(i,j) = 0.0;
}

void PoissonIntegrator::atPoint(apf::Vector3 const& p, double w, double dv)
{
  apf::NewArray<apf::Vector3> gradBF;
  apf::getShapeGrads(e_,p,gradBF);
  apf::Matrix3x3 J;
  apf::getJacobian(me_,p,J);
  double j = apf::getJacobianDeterminant(J,apf::getDimension(me_));
  for (int a=0; a < numLocalEqs_; ++a)
  for (int b=0; b < numLocalEqs_; ++b)
  for (int n=0; n < numDims_; ++n)
    Ke_(a,b) = gradBF[a][n]*gradBF[a][n]*w*j;
}

}

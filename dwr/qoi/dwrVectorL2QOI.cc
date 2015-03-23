/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <apfMesh.h>
#include "dwrVectorL2QOI.h"

namespace dwr {

VectorL2QOI::VectorL2QOI(int o, apf::Field* u) :
  apf::Integrator(o),
  u_(u)
{
  numDims_ = apf::getMesh(u_)->getDimension();
}

VectorL2QOI::~VectorL2QOI()
{
}

void VectorL2QOI::inElement(apf::MeshElement* me)
{
  e_ = apf::createElement(u_,me);
  numNodes_ = apf::countNodes(e_);
  numDofs_ = numDims_ * numNodes_;

  Fe.setSize(numDofs_);
  for (int i=0; i < numDofs_; ++i)
    Fe(i) = 0.0;
}

void VectorL2QOI::outElement()
{
  apf::destroyElement(e_);
}

void VectorL2QOI::atPoint(apf::Vector3 const& p, double w, double dv)
{
  apf::NewArray<double> bf;
  apf::getShapeValues(e_,p,bf);
  apf::Vector3 uVals;
  apf::getVector(e_,p,uVals);
  for (int i=0; i < numNodes_; ++i)
  for (int d=0; d < numDims_; ++d)
    Fe(i*numDims_ + d) += 2.0 * bf[i] * uVals[d] * w * dv;
}

}

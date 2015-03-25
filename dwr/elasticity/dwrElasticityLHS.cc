/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <apfMesh.h>
#include "dwrUtils.h"
#include "dwrGrad.h"
#include "dwrStrain.h"
#include "dwrStress.h"
#include "dwrElasticityLHS.h"

namespace dwr {

ElasticityLHS::ElasticityLHS(int o, apf::Field* u) : 
  apf::Integrator(o),
  primal_(u)
{
  numDims_ = apf::getMesh(primal_)->getDimension();
}

ElasticityLHS::~ElasticityLHS()
{
}

void ElasticityLHS::inElement(apf::MeshElement* me)
{
  e_ = apf::createElement(primal_,me);
  numNodes_ = apf::countNodes(e_);
  numDofs_ = numDims_ * numNodes_;

  u_.setSize(numNodes_);
  residual_.setSize(numDofs_);
  for (int i=0; i < numNodes_; ++i) {
    for (int d=0; d < numDims_; ++d) {
      u_[i][d] = 0.0;
      u_[i][d].diff(i*numDims_ + d, numDofs_);
      residual_[i*numDims_ + d]  = 0.0;
    }
  }

  Ke.setSize(numDofs_,numDofs_);
  for (int i=0; i < numDofs_; ++i)
  for (int j=0; j < numDofs_; ++j)
    Ke(i,j) = 0.0;
}

void ElasticityLHS::outElement()
{
  apf::destroyElement(e_);
}

void ElasticityLHS::atPoint(apf::Vector3 const& p, double w, double dv)
{
  apf::NewArray<apf::Vector3> gradBF;
  apf::getShapeGrads(e_,p,gradBF);

  AD_Matrix3x3 gradU;
  computeVectorGrad(numDims_,numNodes_,u_,gradBF,gradU);

  AD_Matrix3x3 strain;
  computeStrain(gradU,strain);

  AD_Matrix3x3 stress;
  computeStress(E_,nu_,strain,stress);

  for (int i=0; i < numNodes_; ++i)
  for (int c=0; c < numDims_; ++c)
  for (int d=0; d < numDims_; ++d)
    residual_[i*numDims_ + c] += stress[c][d] * gradBF[i][d] * w * dv;

  for (int i=0; i < numDofs_; ++i)
  for (int j=0; j < numDofs_; ++j)
    Ke(i,j) = residual_[i].fastAccessDx(j);
}

}

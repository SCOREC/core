/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <apfMesh.h>
#include <cstring>
#include "dwrElasticity.h"

namespace dwr {

static void zero(AD_Matrix3x3 m)
{
  for (int i=0; i < 3; ++i)
  for (int j=0; j < 3; ++j)
    m[i][j] = 0.0;
}

static void computeGradU(
    int numDims,
    int numNodes,
    apf::DynamicArray<AD_Vector3> const &u,
    apf::NewArray<apf::Vector3> const& gradBF,
    AD_Matrix3x3& gradU)
{
  zero(gradU);
  for (int i=0; i < numNodes; ++i)
  for (int c=0; c < numDims; ++c)
  for (int d=0; d < numDims; ++d)
    gradU[c][d] += u[i][c] * gradBF[i][d];
}

static void computeStrain(AD_Matrix3x3 gradu, AD_Matrix3x3& strain)
{
  zero(strain);
  for (int i=0; i < 3; ++i)
  for (int j=0; j < 3; ++j)
    strain[i][j] = 0.5 * ( gradu[i][j] + gradu[j][i] );
}

static AD trace(AD_Matrix3x3& m)
{
  return m[0][0] + m[1][1] + m[2][2];
}

static void computeStress(AD_Matrix3x3 strain,
    double E, double nu, AD_Matrix3x3& stress)
{
  zero(stress);
  double lm = ( E * nu ) / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu) );
  double mu = E / ( 2.0 * (1.0 + nu) );
  stress[0][0] = 2.0 * mu * strain[0][0] + lm * trace(strain);
  stress[1][1] = 2.0 * mu * strain[1][1] + lm * trace(strain);
  stress[2][2] = 2.0 * mu * strain[2][2] + lm * trace(strain);
  stress[0][1] = 2.0 * mu * strain[0][1];
  stress[1][2] = 2.0 * mu * strain[1][2];
  stress[2][0] = 2.0 * mu * strain[2][0];
  stress[1][0] = stress[0][1];
  stress[2][1] = stress[1][2];
  stress[0][2] = stress[2][0];
}

LinElastInt::LinElastInt(int o, apf::Field* u) : 
  apf::Integrator(o),
  primal_(u)
{
  numDims_ = apf::getMesh(primal_)->getDimension();
}

LinElastInt::~LinElastInt()
{
}

void LinElastInt::inElement(apf::MeshElement* me)
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

void LinElastInt::outElement()
{
  apf::destroyElement(e_);
}

void LinElastInt::atPoint(apf::Vector3 const& p, double w, double dv)
{
  apf::NewArray<apf::Vector3> gradBF;
  apf::getShapeGrads(e_,p,gradBF);

  AD_Matrix3x3 gradU;
  computeGradU(numDims_,numNodes_,u_,gradBF,gradU);

  AD_Matrix3x3 strain;
  computeStrain(gradU,strain);

  AD_Matrix3x3 stress;
  computeStress(strain,E_,nu_,stress);

  for (int i=0; i < numNodes_; ++i)
  for (int c=0; c < numDims_; ++c)
  for (int d=0; d < numDims_; ++d)
    residual_[i*numDims_ + c] += stress[c][d] * gradBF[i][d] * w * dv;

  for (int i=0; i < numDofs_; ++i)
  for (int j=0; j < numDofs_; ++j)
    Ke(i,j) = residual_[i].fastAccessDx(j);

}

}

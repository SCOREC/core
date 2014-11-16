/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrQoI.h"
#include "awrDomainIntegral.h"
#include <apf.h>
#include <apfMesh.h>
#include <apfDynamicVector.h>
#include <Teuchos_ParameterList.hpp>

namespace awr {

/** integrator **/

class DomainIntegrator : public apf::Integrator
{
  public:
    DomainIntegrator(apf::Field* f, int o) :
      apf::Integrator(o),
      f_(f)
    {
      numComponents_ = countComponents(f_);
      assert(numComponents_ == 1);
    }
    void inElement(apf::MeshElement* me)
    {
      me_ = me;
      e_ = apf::createElement(f_,me_);
      numDofs_ = apf::countNodes(e_);
      Fe.setSize(numDofs_);
      for (int a=0; a < numDofs_; ++a)
        Fe(a) = 0.0;
    }
    void outElement()
    {
      apf::destroyElement(e_);
    }
    void atPoint(apf::Vector3 const& p, double w, double dv)
    {
      apf::NewArray<double> bf;
      apf::getShapeValues(e_,p,bf);
      apf::Matrix3x3 J;
      apf::getJacobian(me_,p,J);
      int numDims = apf::getMesh(f_)->getDimension();
      double val = getScalar(e_,p);
      double j = apf::getJacobianDeterminant(J,numDims);
      for (int a=0; a < numDofs_; ++a)
        Fe(a) += bf[a] * val * w * j;
    }
    apf::DynamicVector Fe;
  private:
    apf::Field* f_;
    apf::Element* e_;
    apf::MeshElement* me_;
    int numComponents_;
    int numDofs_;
};

/** qoi **/

DomainIntegral::DomainIntegral(
    ParameterList& p, apf::Mesh* m, apf::Field* f) :
  QoI(p,m,f)
{
}

DomainIntegral::~DomainIntegral()
{
  delete integrator_;
}

void DomainIntegral::validateQoIList()
{
}

void DomainIntegral::createIntegrator()
{
  integrator_ = new DomainIntegrator(primal_,integrationOrder_);
}

void DomainIntegral::processFe(
    apf::MeshEntity* e,
    apf::DynamicVector& Fe)
{
  apf::MeshElement* me = apf::createMeshElement(mesh_,e);
  integrator_->process(me);
  Fe = integrator_->Fe;
  apf::destroyMeshElement(me);
}

}

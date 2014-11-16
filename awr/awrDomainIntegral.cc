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
    }
    void inElement(apf::MeshElement* me)
    {
    }
    void outElement()
    {
    }
    void atPoint(apf::Vector3 const& p, double w, double dv)
    {
    }
    apf::DynamicVector Fe;
  private:
    apf::Field* f_;
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

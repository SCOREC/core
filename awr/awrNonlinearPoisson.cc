/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrProblem.h"
#include "awrNonlinearPoisson.h"
#include <apf.h>
#include <apfMesh.h>
#include <apfDynamicMatrix.h>
#include <Teuchos_ParameterList.hpp>

namespace awr {

/** integrator **/

class NonlinearPoissonIntegrator : public apf::Integrator
{
  public:
    NonlinearPoissonIntegrator(apf::Field* f, int o) :
      apf::Integrator(o),
      u_(f)
    {
      numDims_ = apf::getMesh(u_)->getDimension();
      int numComponents = apf::countComponents(u_);
      assert(numComponents == 1);
    }
    void inElement(apf::MeshElement* me)
    {
      me_ = me;
      e_ = apf::createElement(u_,me_);
      numNodes_ = apf::countNodes(e_);
      Ke.setSize(numNodes_,numNodes_);
      for (int a=0; a < numNodes_; ++a)
      for (int b=0; b < numNodes_; ++b)
        Ke(a,b) = 0.0;
    }
    void outElement()
    {
      apf::destroyElement(e_);
    }
    void atPoint(apf::Vector3 const& p, double w, double dv)
    {
      apf::NewArray<double> bf;
      apf::getShapeValues(e_,p,bf);
      apf::NewArray<apf::Vector3> gradBF;
      apf::getShapeGrads(e_,p,gradBF);
      double u = getScalar(e_,p);
      apf::Vector3 gradU;
      getGrad(e_,p,gradU);
      for (int a=0; a < numNodes_; ++a)
      for (int b=0; b < numNodes_; ++b)
      for (int i=0; i < numDims_; ++i)
        Ke(a,b) += (2.0 * u * gradU[i] * bf[a] * gradBF[b][i] +
            (1.0 + u*u) * gradBF[a][i] * gradBF[a][i] ) * w * dv;
    }
    apf::DynamicMatrix Ke;
  private:
    apf::Field* u_;
    apf::Element* e_;
    apf::MeshElement* me_;
    int numDims_;
    int numNodes_;
};

/** problem **/

NonlinearPoissonProblem::
NonlinearPoissonProblem(ParameterList& p, apf::Mesh* m) :
  Problem(p,m)
{
}

NonlinearPoissonProblem::~NonlinearPoissonProblem()
{
  delete integrator_;
}

void rejectNonlinearPoissonInput(const char* msg)
{
  fprintf(stderr,"AWR Nonlinear Poisson problem input error\n");
  fprintf(stderr,"%s\n",msg);
  abort();
}

void NonlinearPoissonProblem::validateProblemList()
{
  if (! problemList_.isParameter("Primal Solution Field Name"))
    rejectNonlinearPoissonInput("*Primal Solution Field Name* not set");
}

void NonlinearPoissonProblem::setPrimalField()
{
  std::string name = problemList_.get("Primal Solution Field Name","");
  primal_ = mesh_->findField(name.c_str());
  if (primal_ == NULL)
    rejectNonlinearPoissonInput("Invalid primal solution field name");
}

void NonlinearPoissonProblem::createIntegrator()
{
  integrator_ = new NonlinearPoissonIntegrator(primal_,integrationOrder_);
}

void NonlinearPoissonProblem::processKe(
    apf::MeshEntity* e,
    apf::DynamicMatrix& Ke)
{
  apf::MeshElement* me = apf::createMeshElement(mesh_,e);
  integrator_->process(me);
  Ke = integrator_->Ke;
  apf::destroyMeshElement(me);
}

}

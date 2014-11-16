/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrProblem.h"
#include "awrPoisson.h"
#include <apf.h>
#include <apfMesh.h>
#include <apfDynamicMatrix.h>
#include <Teuchos_ParameterList.hpp>

namespace awr {

/** integrator **/

class PoissonIntegrator : public apf::Integrator
{
  public:
    PoissonIntegrator(apf::Field* f, int o) :
      apf::Integrator(o),
      f_(f)
    {
      numDims_ = apf::getMesh(f_)->getDimension();
      int numComponents = apf::countComponents(f_);
      assert(numComponents == 1);
    }
    void inElement(apf::MeshElement* me)
    {
      me_ = me;
      e_ = apf::createElement(f_,me_);
      numNodes = apf::countNodes(e_);
      Ke.setSize(numNodes,numNodes);
      for (int a=0; a < numNodes; ++a)
      for (int b=0; b < numNodes; ++b)
        Ke(a,b) = 0.0;
    }
    void outElement()
    {
      apf::destroyElement(e_);
    }
    void atPoint(apf::Vector3 const& p, double w, double dv)
    {
      apf::NewArray<apf::Vector3> gradBF;
      apf::getShapeGrads(e_,p,gradBF);
      apf::Matrix3x3 J;
      apf::getJacobian(me_,p,J);
      double j = apf::getJacobianDeterminant(J,numDims_);
      for (int a=0; a < numNodes; ++a)
      for (int b=0; b < numNodes; ++b)
      for (int i=0; i < numDims_; ++i)
        Ke(a,b) += gradBF[a][i] * gradBF[b][i] * w * j;
    }
    apf::DynamicMatrix Ke;
    int numNodes;
  private:
    apf::Field* f_;
    apf::Element* e_;
    apf::MeshElement* me_;
    int numDims_;
};

/** problem **/

PoissonProblem::PoissonProblem(ParameterList& p, apf::Mesh* m) :
  Problem(p,m)
{
}

PoissonProblem::~PoissonProblem()
{
  delete integrator_;
}

void rejectPoissonInput(const char* msg)
{
  fprintf(stderr,"AWR Poisson problem input error\n");
  fprintf(stderr,"%s\n",msg);
  abort();
}

void PoissonProblem::validateProblemList()
{
  if (! problemList_.isParameter("Primal Solution Field Name"))
    rejectPoissonInput("*Primal Solution Field Name* not set");
}

void PoissonProblem::setPrimalField()
{
  std::string name = problemList_.get("Primal Solution Field Name","");
  primal_ = mesh_->findField(name.c_str());
  if (primal_ == NULL)
    rejectPoissonInput("Invalid primal solution field name");
}

void PoissonProblem::createIntegrator()
{
  integrator_ = new PoissonIntegrator(primal_,integrationOrder_);
}

void PoissonProblem::processKe(
    apf::MeshEntity* e,
    apf::DynamicMatrix& Ke)
{
  apf::MeshElement* me = apf::createMeshElement(mesh_,e);
  integrator_->process(me);
  Ke = integrator_->Ke;
  apf::destroyMeshElement(me);
}

}

/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "../awrUtils.h"
#include "awrPoissonProblem.h"
#include <apfMesh.h>

namespace awr {

PoissonProblem::PoissonProblem(apf::Mesh* m,
    const Teuchos::RCP<Teuchos::ParameterList>& p):
  Problem(m,p)
{
  validateParameters();
  numEquations_ = 1;
}

void checkFieldExists(std::string name, apf::Field* f)
{
  if (f == NULL)
    fail((name + " does not exist").c_str());
}

void PoissonProblem::validateParameters()
{
  Teuchos::ParameterList& problemList =
    params_->sublist("Adjoint Problem");
  integrationOrder_ = problemList.get("Integration Order",2);
  std::string name = problemList.get("Primal Solution Field Name","");
  primalSolution_ = mesh_->findField(name.c_str());
  checkFieldExists(name,primalSolution_);
}

}

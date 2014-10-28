/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "../awrUtils.h"
#include "awrProblemFactory.h"
#include "awrPoissonProblem.h"

namespace awr {

ProblemFactory::ProblemFactory(apf::Mesh* m,
    const Teuchos::RCP<Teuchos::ParameterList>& p):
  mesh_(m),
  params_(p)
{
}

Teuchos::RCP<Problem> ProblemFactory::create()
{
  Teuchos::RCP<Problem> strategy;
  Teuchos::ParameterList problem_sublist =
    params_->sublist("Adjoint Problem",false);
  std::string method = problem_sublist.get("Name","");
  if (method == "Poisson")
    strategy = Teuchos::rcp(new PoissonProblem(mesh_,params_));
  else
    fail((method + " is not a valid problem").c_str());
  return strategy;
}

}

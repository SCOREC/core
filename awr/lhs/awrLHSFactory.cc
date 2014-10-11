/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrLHSFactory.h"
#include "awrPoissonLHS.h"
#include "awrNonlinearPoissonLHS.h"

namespace awr {

LHSFactory::LHSFactory(apf::Mesh* m, const Teuchos::ParameterList& p) :
  mesh_(m),
  params_(p)
{
}

Teuchos::RCP<LHS> LHSFactory::create()
{
  Teuchos::RCP<LHS> strategy;
  std::string method = params_.get("Adjoint Problem Name","");
  if (method == "Poisson")
    strategy = Teuchos::rcp(new PoissonLHS(mesh_,params_));
  else if (method == "Nonlinear Poisson")
    strategy = Teuchos::rcp(new NonlinearPoissonLHS(mesh_,params_));
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, Teuchos::Exceptions::InvalidParameter,
        "AWR: LHS Factory: Adjoint Problem Name "
        << method << " does not exist");
  return strategy;
}

}

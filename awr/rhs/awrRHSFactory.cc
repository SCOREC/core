/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrRHSFactory.h"
#include "awrPoissonRHS.h"
#include "awrNonlinearPoissonRHS.h"

namespace awr {

RHSFactory::RHSFactory(apf::Mesh* m, const Teuchos::ParameterList& p) :
  mesh_(m),
  params_(p)
{
}

Teuchos::RCP<RHS> RHSFactory::create()
{
  Teuchos::RCP<RHS> strategy;
  std::string method = params_.get("Adjoint Problem Name","");
  if (method == "Poisson")
    strategy = Teuchos::rcp(new PoissonRHS(mesh_,params_));
  else if (method == "Nonlinear Poisson")
    strategy = Teuchos::rcp(new NonlinearPoissonRHS(mesh_,params_));
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, Teuchos::Exceptions::InvalidParameter,
        "AWR: RHS Factory: Adjoint Problem Name "
        << method << " does not exist");
  return strategy;
}

}

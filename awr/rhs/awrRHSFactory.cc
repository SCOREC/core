/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apf.h"
#include "awrRHSFactory.h"
#include "awrNonlinearPoissonRHS.h"

namespace awr {

RHSFactory::RHSFactory(const Teuchos::ParameterList& p) :
  params_(p)
{
}

Teuchos::RCP<RHS> RHSFactory::create()
{
  Teuchos::RCP<RHS> strategy;
  std::string method = params_.get("Adjoint Problem Name","");
  if (method == "Nonlinear Poisson")
    strategy = Teuchos::rcp(new NonlinearPoissonRHS(params_));
  else
    apf::fail("AWR: unknown adjoint problem name");
  return strategy;
}

}

/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrQOIFactory.h"
#include "awrSolutionSquaredQOI.h"

namespace awr {

QOIFactory::QOIFactory(apf::Mesh* m, const Teuchos::ParameterList& p) :
  mesh_(m),
  params_(p)
{
}

Teuchos::RCP<QOI> QOIFactory::create()
{
  Teuchos::RCP<QOI> strategy;
  std::string method = params_.get("Quantity of Interest Name","");
  if (method == "Solution Squared")
    strategy = Teuchos::rcp(new SolutionSquaredQOI(mesh_,params_));
  else
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, Teuchos::Exceptions::InvalidParameter,
        "AWR: QOI Factory: Quantity of Interest Name "
        << method << " does not exist");
  return strategy;
}

}

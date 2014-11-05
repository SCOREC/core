/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "../awrUtils.h"
#include "awrProblem.h"

namespace awr {

void validate(std::string name, bool valid)
{
  if (!valid)
    fail((name + " is not a sublist").c_str());
}

void Problem::validateSublists()
{
  bool problemDefined = params_->isSublist("Adjoint Problem");
  validate("Adjoint Problem",problemDefined);
  bool qoiDefined = params_->isSublist("Quantity Of Interest");
  validate("Quantity of Interest",qoiDefined);
  bool bcDefined = params_->isSublist("Boundary Conditions");
  validate("Boundary Conditions",bcDefined);
  bool solverDefined = params_->isSublist("Solver");
  validate("Solver",solverDefined);
}

Problem::Problem(apf::Mesh* m,
    const Teuchos::RCP<Teuchos::ParameterList>& p):
  mesh_(m),
  params_(p)
{
  validateSublists();
}

}

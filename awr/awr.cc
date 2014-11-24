/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awr.h"
#include "awrProblem.h"
#include <Teuchos_ParameterList.hpp>

namespace awr {

void solveAdjoint(ParameterList& p, apf::Mesh* m)
{
  Problem* problem = createProblem(p,m);
  problem->setup();
  problem->assemble();
  problem->solve();
  delete problem;
}

}

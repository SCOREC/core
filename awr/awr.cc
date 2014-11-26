/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awr.h"
#include "awrProblem.h"
#include <PCU.h>
#include <Teuchos_ParameterList.hpp>

namespace awr {

void solveAdjoint(ParameterList& p, apf::Mesh* m)
{
  double t0 = MPI_Wtime();
  Problem* problem = createProblem(p,m);
  problem->setup();
  problem->assemble();
  problem->solve();
  delete problem;
  double t1 = MPI_Wtime();
  print("solveAdjoint took %f seconds total\n",t1-t0);
}

}

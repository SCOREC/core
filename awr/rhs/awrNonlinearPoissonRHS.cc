/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrNonlinearPoissonRHS.h"

namespace awr {

NonlinearPoissonRHS::NonlinearPoissonRHS(const Teuchos::ParameterList& p) :
  RHS(p)
{
}

void NonlinearPoissonRHS::evaluateElementRHS()
{
}

}

/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrPoissonRHS.h"

namespace awr {

/*****************************************************************************/
PoissonRHS::
PoissonRHS(const Teuchos::ParameterList& p) :
  RHS(p)
{
}

/*****************************************************************************/
void
PoissonRHS::
evaluateElementRHS(apf::MeshEntity* element,
                   apf::Field* primal_solution,
                   int integration_order,
                   apf::DynamicMatrix& k)
{
}

/*****************************************************************************/
}

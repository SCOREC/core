/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrNonlinearPoissonRHS.h"

namespace awr {

/*****************************************************************************/
NonlinearPoissonRHS::
NonlinearPoissonRHS(apf::Mesh* m, const Teuchos::ParameterList& p) :
  RHS(m,p)
{
}

/*****************************************************************************/
void
NonlinearPoissonRHS::
evaluateElementRHS(apf::MeshEntity* element,
                   int integration_order,
                   apf::DynamicMatrix& k)
{
}

/*****************************************************************************/
}

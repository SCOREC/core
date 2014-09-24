/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrPoissonRHS.h"
#include "Teuchos_TestForException.hpp"

namespace awr {

/*****************************************************************************/
PoissonRHS::
PoissonRHS(apf::Mesh* m, const Teuchos::ParameterList& p) :
  RHS(m,p)
{
}

/*****************************************************************************/
void
PoissonRHS::
evaluateElementRHS(apf::MeshEntity* element,
                   int integration_order,
                   apf::DynamicMatrix& k)
{
}

/*****************************************************************************/
void
PoissonRHS::
validateParameters()
{
  std::string n = params_.get("Primal Solution Field Name","");
}

/*****************************************************************************/
}

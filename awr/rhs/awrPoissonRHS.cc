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
  validateParameters();
}

/*****************************************************************************/
void
PoissonRHS::
evaluateElementRHS(apf::MeshEntity* element,
                   apf::DynamicMatrix& k)
{
}

/*****************************************************************************/
void
PoissonRHS::
validateParameters()
{
  integration_order_ = params_.get("Integration Order",2); 
  std::string name = params_.get("Primal Solution Field Name","");
  sol_ = mesh_->findField(name.c_str());
  if (sol_ == NULL)
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, Teuchos::Exceptions::InvalidParameter,
        "AWR: Poisson RHS: solution field with name \""
        << name << "\" does not exist\n");
                             
}

/*****************************************************************************/
}

/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrNonlinearPoissonLHS.h"

namespace awr {

/*****************************************************************************/
NonlinearPoissonLHS::
NonlinearPoissonLHS(apf::Mesh* m, const Teuchos::ParameterList& p) :
  LHS(m,p)
{
}

/*****************************************************************************/
void
NonlinearPoissonLHS::
evaluateElementLHS(apf::MeshEntity* element,
                   apf::DynamicMatrix& k)
{
}

/*****************************************************************************/
}

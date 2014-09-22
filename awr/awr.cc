/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awr.h"
#include "apfMesh.h"
#include "apfField.h"
#include "apfShape.h"
#include "rhs/awrRHS.h"
#include "rhs/awrRHSFactory.h"

namespace awr {

apf::Field* enrichSolution(apf::Field* sol, const char* name_e)
{
  apf::Mesh* m = sol->getMesh();
  std::string name = apf::getName(sol);
  int type = apf::getValueType(sol);
  apf::FieldShape* shape = sol->getShape();
  assert (shape == apf::getLagrange(1));
  int order = shape->getOrder();
  int order_e = order+1;
  apf::FieldShape* shape_e = apf::getLagrange(order_e);
  apf::Field* sol_e =
    apf::createField(m,name_e,type,shape_e);
  sol_e->project(sol);
  apf::destroyField(sol);
  apf::changeMeshShape(m,shape_e,/*project=*/true);
  return sol_e;
}

apf::Field* solveAdjointProblem(const Teuchos::ParameterList& params)
{
  RHSFactory rhsFactory(params);
  Teuchos::RCP<RHS> rhs = rhsFactory.create();
}

}

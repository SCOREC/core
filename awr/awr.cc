/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awr.h"
#include "awrProblem.h"
#include <PCU.h>
#include <apf.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <Teuchos_ParameterList.hpp>

namespace awr {


apf::Field* enrichSolution(apf::Field* sol, const char* name_e)
{
  apf::Mesh* m = apf::getMesh(sol);
  std::string name = apf::getName(sol);
  int type = apf::getValueType(sol);
  apf::FieldShape* shape = apf::getShape(sol);
  assert (shape == apf::getLagrange(1));
  int order = shape->getOrder();
  int order_e = order+1;
  apf::FieldShape* shape_e = apf::getLagrange(order_e);
  apf::Field* sol_e =
    apf::createField(m,name_e,type,shape_e);
  apf::projectField(sol_e,sol);
  apf::destroyField(sol);
  apf::changeMeshShape(m,shape_e,/*project=*/true);
  return sol_e;
}

void solveAdjoint(ParameterList& p, apf::Mesh* m)
{
  double t0 = PCU_Time();
  Problem* problem = createProblem(p,m);
  problem->setup();
  problem->assemble();
  problem->solve();
  delete problem;
  double t1 = PCU_Time();
  print("solveAdjoint took %f seconds total\n",t1-t0);
}

}

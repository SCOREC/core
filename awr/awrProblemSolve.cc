/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrProblem.h"
#include "awrLinearSystem.h"
#include <apfNumbering.h>
#include <iostream>

namespace awr {

void attachSolution(
    apf::Mesh* m,
    apf::GlobalNumbering* gn,
    apf::Field* f,
    double* sol)
{
  apf::DynamicArray<apf::Node> nodes;
  getNodes(gn,nodes);
  int nc = apf::countComponents(f);
  double v[nc];
  for (int i=0; i < nodes.getSize(); ++i)
  {
    for (int c=0; c < nc; ++c)
    {
      v[c] = sol[i+c];
    }
    apf::setComponents(f,nodes[i].entity,nodes[i].node,v);
  }
}

void Problem::solve()
{
  ls_->solve();
  double* sol = ls_->getSolution();
  attachSolution(mesh_,globalNumbering_,adjoint_,sol);
}

}

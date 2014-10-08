/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrRHS.h"
#include <apfMesh.h>

namespace awr {

RHS::RHS(apf::Mesh* m, const Teuchos::ParameterList& p) :
  mesh_(m),
  params_(p)
{
}

void RHS::assemble()
{
  apf::MeshIterator* elems = mesh_->begin(mesh_->getDimension());
  apf::MeshEntity* e = mesh_->iterate(elems);
  mesh_->end(elems);
  apf::DynamicMatrix k;
  evaluateElementRHS(e,k);
  for (int i=0; i < k.getRows(); ++i)
  {
    for (int j=0; j < k.getColumns(); ++j)
      std::cout << k(i,j) << " ";
    std::cout << std::endl;
  }
}

}

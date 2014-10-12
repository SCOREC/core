/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrQOI.h"
#include <apfMesh.h>

namespace awr {

QOI::QOI(apf::Mesh* m, const Teuchos::ParameterList& p) :
  mesh_(m),
  params_(p)
{
}

void QOI::assemble()
{
  apf::MeshIterator* elems = mesh_->begin(mesh_->getDimension());
  apf::MeshEntity* e = mesh_->iterate(elems);
  mesh_->end(elems);
  apf::DynamicVector k;
  evaluateElementQOI(e,k);
}

}

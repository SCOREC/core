/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrRHS.h"

namespace awr {

RHS::RHS(apf::Mesh* m, const Teuchos::ParameterList& p) :
  mesh_(m),
  params_(p)
{
}

void RHS::assemble()
{
}

}

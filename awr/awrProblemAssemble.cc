/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "awrProblem.h"
#include "awrLinearSystem.h"
#include <apf.h>
#include <apfMesh.h>
#include <apfNumbering.h>
#include <apfDynamicMatrix.h>
#include <Teuchos_ParameterList.hpp>

namespace awr {

void Problem::processBC()
{
}

void Problem::assemble()
{
  /* pure virtual method */
  createIntegrator();
  processBC();
}

}

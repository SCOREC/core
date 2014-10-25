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

}

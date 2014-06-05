/*
 * Copyright (C) 2014 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfZoltan.h"
#include <apf.h>

namespace apf {

Splitter* makeZoltanSplitter(Mesh* mesh, int method, int approach,
    bool debug, bool sync)
{
  fail("apf_zoltan compiled empty !");
  return 0;
}

Balancer* makeZoltanBalancer(Mesh* mesh, int method, int approach, bool debug)
{
  fail("apf_zoltan compiled empty !");
  return 0;
}

}

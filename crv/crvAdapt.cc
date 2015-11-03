/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <apf.h>
#include <PCU.h>
#include "crvAdapt.h"
#include <maBalance.h>
#include <maCoarsen.h>
#include <maSnap.h>

namespace crv {

void uniformRefine(ma::Mesh* m)
{
  ma::Input* in = ma::configureUniformRefine(m, 1);
  if (in->shouldSnap) {
    in->shouldSnap = false;
  }
  in->shouldFixShape = false;
  in->shapeHandler = getShapeHandler(m);
  crv::adapt(in);
}

void adapt(ma::Input* in)
{
  ma::print("version 2.0 !");
  double t0 = PCU_Time();
  ma::validateInput(in);
  ma::Adapt* a = new ma::Adapt(in);
  ma::preBalance(a);
  for (int i=0; i < in->maximumIterations; ++i)
  {
    ma::print("iteration %d",i);
    ma::coarsen(a);
    ma::midBalance(a);
    crv::refine(a);
  }
  ma::snap(a);
  ma::postBalance(a);
  delete a;
  delete in;
  double t1 = PCU_Time();
  ma::print("mesh adapted in %f seconds",t1-t0);
  apf::printStats(a->mesh);
}

}

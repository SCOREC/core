/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crvAdapt.h"
#include "maTemplates.h"
#include <PCU.h>
#include <cassert>

namespace crv {



bool refine(ma::Adapt* a)
{
  double t0 = PCU_Time();
  --(a->refinesLeft);
  long count = ma::markEdgesToSplit(a);
  if ( ! count) {
    return false;
  }
  assert(ma::checkFlagConsistency(a,1,ma::SPLIT));
  ma::Refine* r = a->refine;
  ma::resetCollection(r);
  ma::collectForTransfer(r);
  ma::collectForMatching(r);
  ma::addAllMarkedEdges(r);
  ma::splitElements(r);
  ma::processNewElements(r);
  ma::destroySplitElements(r);
  crv::repositionInterior(r);
  ma::forgetNewEntities(r);
  double t1 = PCU_Time();
  ma::print("refined %li edges in %f seconds",count,t1-t0);
  return true;
}

}

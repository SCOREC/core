#include "phAdapt.h"
#include "chef.h"
#include "ph.h"
#include <ma.h>
#include <PCU.h>
#include <sam.h>

#include <cstdio>
#include <cassert>

namespace ph {

void setupMatching(ma::Input& in) {
  if (!PCU_Comm_Self())
    printf("Matched mesh: disabling"
           " snapping, and shape correction,\n");
  in.shouldSnap = false;
  in.shouldFixShape = false;
}

static void runFromErrorThreshold(Input&, apf::Mesh2* m)
{
  const unsigned idx = 5;
  const double errLimit = 1e-6;
  const double factor = 0.5;
  apf::Field* szFld = sam::errorThreshold(m,"errors",idx,errLimit,factor);
  assert(szFld);
  chef::adapt(m, szFld);
  apf::destroyField(szFld);
}

static void runFromGivenSize(Input&, apf::Mesh2* m)
{
  const unsigned idx = 5;
  apf::Field* szFld = sam::specifiedIso(m,"errors",idx);
  assert(szFld);
  chef::adapt(m, szFld);
  apf::destroyField(szFld);
}

void tetrahedronize(Input&, apf::Mesh2* m)
{
  ma::Input* ma_in = ma::configureIdentity(m);
  ma_in->shouldRunPreParma = true;
  ma_in->shouldTurnLayerToTets = true;
  ma::adapt(ma_in);
  m->verify();
}

void adapt(Input& in, apf::Mesh2* m)
{
  typedef void (*Strategy)(Input&, apf::Mesh2*);
  static Strategy const table[PH_STRATEGIES] =
  {0//0
  ,runFromGivenSize//1
  ,runFromErrorThreshold//2
  ,0//3
  ,0//4
  ,0//5
  ,0//6
  ,chef::uniformRefinement //7
  ,0//8
  };
  table[in.adaptStrategy](in, m);
  m->verify();
}

}

namespace chef {
  void adapt(apf::Mesh2* m, apf::Field* szFld) {
    ma::Input* ma_in = ma::configure(m, szFld);
    ma_in->shouldRunPreZoltan = true;
    if (m->hasMatching())
      ph::setupMatching(*ma_in);
    ma::adapt(ma_in);
  }

  void uniformRefinement(ph::Input& in, apf::Mesh2* m)
  {
    ma::Input* ma_in = ma::configureMatching(m, in.recursiveUR);
    ma_in->shouldRefineLayer = true;
    ma_in->splitAllLayerEdges = in.splitAllLayerEdges;
    if (in.snap) {
      if (!ma_in->shouldSnap)
        ph::fail("adapt.inp requests snapping but model doesn't support it\n");
    } else
      ma_in->shouldSnap = false;
    ma::adapt(ma_in);
  }
}

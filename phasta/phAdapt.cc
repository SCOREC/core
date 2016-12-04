#include "phAdapt.h"
#include "chef.h"
#include "ph.h"
#include <ma.h>
#include <PCU.h>
#include <sam.h>

#include <cstdio>
#include <cassert>

namespace ph {

void setupPreBalance(Input& in, ma::Input* ma_in) {
  if ( in.preAdaptBalanceMethod == "parma" ) {
    ma_in->shouldRunPreParma = true;
  } else if( in.preAdaptBalanceMethod == "graph" ) {
    ma_in->shouldRunPreZoltan = true;
  } else if( in.preAdaptBalanceMethod == "zrib" ) {
    ma_in->shouldRunPreZoltanRib = true;
  } else if ( in.preAdaptBalanceMethod == "none" ) {
    ma_in->shouldRunPreZoltan = false;
    ma_in->shouldRunPreZoltanRib = false;
    ma_in->shouldRunPreParma = false;
  } else {
    if (!PCU_Comm_Self())
      fprintf(stderr,
          "warning: ignoring unknown value of preAdaptBalanceMethod %s\n",
          in.preAdaptBalanceMethod.c_str());
  }
}

void setupMatching(ma::Input& in) {
  if (!PCU_Comm_Self())
    printf("Matched mesh: disabling"
           " snapping, and shape correction,\n");
  in.shouldSnap = false;
  in.shouldFixShape = false;
}

static void runFromErrorThreshold(Input& in, apf::Mesh2* m)
{
  const char* fieldname = in.adaptErrorFieldName.c_str();
  const unsigned idx = in.adaptErrorFieldIndex;
  const double errLimit = in.adaptErrorThreshold;
  const double factor = 0.5;
  apf::Field* szFld = sam::errorThreshold(m,fieldname,idx,errLimit,factor);
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

void tetrahedronize(Input& in, apf::Mesh2* m)
{
  ma::Input* ma_in = ma::configureIdentity(m);
  setupPreBalance(in, ma_in);
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

  void adapt(apf::Mesh2* m, apf::Field* szFld, ph::Input& in) {
    ma::Input* ma_in = ma::configure(m, szFld);
    ma_in->shouldRunPreZoltan = true;
    ma_in->shouldTransferParametric = in.transferParametric;
    ma_in->shouldRunMidParma = true; 
    ma_in->shouldRunPostParma = true; 
    ma_in->shouldSnap = in.snap;
    ma_in->maximumIterations = in.maxAdaptIterations;
    if (m->hasMatching())
      ph::setupMatching(*ma_in);
    ma::adapt(ma_in);
  }

  void uniformRefinement(ph::Input& in, apf::Mesh2* m)
  {
    ma::Input* ma_in = ma::configureMatching(m, in.recursiveUR);
    setupPreBalance(in, ma_in);
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

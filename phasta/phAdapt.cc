#include "phAdapt.h"
#include <ma.h>

namespace ph {

static void runUniformRefinement(Input& in, apf::Mesh2* m)
{
  ma::Input* ma_in = ma::configureMatching(m, in.recursiveUR);
  ma_in->shouldRefineLayer = true;
  ma::adapt(ma_in);
}

void tetrahedronize(Input&, apf::Mesh2* m)
{
  ma::Input* ma_in = ma::configureIdentity(m);
  ma_in->shouldRunPreParma = true;
  ma_in->shouldCleanupLayer = true;
  ma_in->shouldTurnLayerToTets = true;
  ma::adapt(ma_in);
  m->verify();
}

void adapt(Input& in, apf::Mesh2* m)
{
  typedef void (*Strategy)(Input&, apf::Mesh2*);
  static Strategy const table[PH_STRATEGIES] = 
  {0//0
  ,0//1
  ,0//2
  ,0//3
  ,0//4
  ,0//5
  ,0//6
  ,runUniformRefinement //7
  ,0//8
  };
  table[in.adaptStrategy](in, m);
  m->verify();
}

}

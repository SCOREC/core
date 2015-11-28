#include "phAdapt.h"
#include "ph.h"
#include <ma.h>
#include <PCU.h>
#include <sam.h>

#include <cstdio>
#include <cassert>

namespace ph {

static void runUniformRefinement(Input& in, apf::Mesh2* m)
{
  ma::Input* ma_in = ma::configureMatching(m, in.recursiveUR);
  ma_in->shouldRefineLayer = true;
  ma_in->splitAllLayerEdges = in.splitAllLayerEdges;
  if (in.snap) {
    if (!ma_in->shouldSnap)
      fail("adapt.inp requests snapping but model doesn't support it\n");
  } else
    ma_in->shouldSnap = false;
  ma::adapt(ma_in);
}

class ReturnErrorSize : public ma::IsotropicFunction
{
  public:
    ReturnErrorSize(apf::Mesh* m)
    {
      const unsigned idx = 5;
      const double errLimit = 1e-6;
      const double factor = 0.5;
      szFld = sam::specifiedIso(m,"errors",idx,errLimit,factor);
      assert(szFld);
    }
    ~ReturnErrorSize()
    {
      apf::destroyField(szFld);
    }
    virtual double getValue(ma::Entity* vert)
    {
      return apf::getScalar(szFld,vert,0);
    }
  private:
    apf::Field* szFld;
};

static void runFromErrorSize(Input& in, apf::Mesh2* m)
{
  ReturnErrorSize sf(m);
  ma::Input* ma_in = ma::configure(m, &sf);
  ma_in->shouldRunPreZoltan = true;
//  ma_in->shouldSnap = false;  //FIXME - is this needed?
  if (m->hasMatching()) {
    if (in.snap)
      fail("adapt.inp requests snapping but mesh is periodic\n");
    if (!PCU_Comm_Self())
      printf("Matched mesh: disabling coarsening, snapping, and shape correction,\n"
             "  synchronizing \"errors\" field (source of size)\n");
    ma_in->shouldCoarsen = false;
    ma_in->shouldSnap = false;
    ma_in->shouldFixShape = false;
    apf::synchronize(m->findField("errors"));
  }
  ma::adapt(ma_in);
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
  ,runFromErrorSize//1
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

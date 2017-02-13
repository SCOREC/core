#include "phAdapt.h"
#include "chef.h"
#include "ph.h"
#include <ma.h>
#include <PCU.h>
#include <sam.h>
#include <parma.h>

#include <cstdio>
#include <cassert>

namespace ph {

struct AdaptCallback : public Parma_GroupCode
{
  apf::Mesh2* mesh;
  ma::Input* in;
  void run(int) {
    ma::adapt(in);
  }
};

static double getAveragePartDensity(apf::Mesh* m) {
  double nElements = m->count(m->getDimension());
  nElements = PCU_Add_Double(nElements);
  return nElements / PCU_Comm_Peers();
}

static int getShrinkFactor(apf::Mesh* m, double minPartDensity) {
  double partDensity = getAveragePartDensity(m);
  int factor = 1;
  while (partDensity < minPartDensity) {
    if (factor >= PCU_Comm_Peers())
      break;
    factor *= 2;
    partDensity *= 2;
  }
  assert(PCU_Comm_Peers() % factor == 0);
  return factor;
}

static void warnAboutShrinking(int factor) {
  int nprocs = PCU_Comm_Peers() / factor;
  if (!PCU_Comm_Self()) {
    fprintf(stderr,"sensing mesh is spread too thin: "
                   "adapting with %d procs\n", nprocs);
  }
}

void adaptShrunken(apf::Mesh2* m, double minPartDensity,
                   Parma_GroupCode& callback) {
  int factor = getShrinkFactor(m, minPartDensity);
  if (!PCU_Comm_Self())
    fprintf(stderr,"adaptShrunken factor computed as %d\n", factor);
  if (factor == 1) {
    callback.run(0);
  } else {
    warnAboutShrinking(factor);
    Parma_ShrinkPartition(m, factor, callback);
  }
}

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

void setupMidBalance(Input& in, ma::Input* ma_in) {
  if ( in.midAdaptBalanceMethod == "parma" ) {
    ma_in->shouldRunMidParma = true;
  } else if( in.midAdaptBalanceMethod == "graph" ) {
    ma_in->shouldRunMidZoltan = true;
  } else if ( in.midAdaptBalanceMethod == "none" ) {
    ma_in->shouldRunMidZoltan = false;
    ma_in->shouldRunMidParma = false;
  } else {
    if (!PCU_Comm_Self())
      fprintf(stderr,
          "warning: ignoring unknown value of midAdaptBalanceMethod %s\n",
          in.midAdaptBalanceMethod.c_str());
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

static void runFromGivenSize(Input& in, apf::Mesh2* m)
{
  const unsigned idx = 5;
  apf::Field* szFld = sam::specifiedIso(m,"errors",idx);
  assert(szFld);
  chef::adapt(m, szFld,in);
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
    ph::AdaptCallback acb;
    acb.mesh = m;
    acb.in = ma_in;
    adaptShrunken(m,10000,acb);
  }

  void adapt(apf::Mesh2* m, apf::Field* szFld, ph::Input& in) {
    ma::Input* ma_in = ma::configure(m, szFld);
    //chef defaults
    ma_in->shouldRunPreZoltan = true;
    ma_in->shouldRunMidParma = true;
    ma_in->shouldRunPostParma = true;
    //override with user inputs if specified
    setupPreBalance(in, ma_in);
    setupMidBalance(in, ma_in);
    ma_in->shouldTransferParametric = in.transferParametric;
    ma_in->shouldSnap = in.snap;
    ma_in->maximumIterations = in.maxAdaptIterations;
    /*
      validQuality sets which elements will be accepted during mesh
      modification. If no boundary layers, you might bring this high (e.g,
      O(1e-2).  This would be bad for boundary layers though since their
      high aspect ratio routinely produces quality measures in the e-6 to e-7
      range so, when there are layers, this needs to be O(1e-8).
    */
    ma_in->validQuality = in.validQuality;
    if (m->hasMatching())
      ph::setupMatching(*ma_in);
    ph::AdaptCallback acb;
    acb.mesh = m;
    acb.in = ma_in;
    adaptShrunken(m,in.adaptShrinkLimit,acb);
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

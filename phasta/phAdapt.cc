#include "phAdapt.h"
#include "chef.h"
#include "ph.h"
#include <ma.h>
#include <PCU.h>
#include <sam.h>
#include <parma.h>

#include <cstdio>
#include <pcu_util.h>

namespace ph {

void setupBalance(const char* key, std::string& method,
    bool& parmaBal, bool& zoltanBal, bool& zoltanRibBal) {
  if ( method == "parma" ) {
    parmaBal = true;
    zoltanBal = false;
    zoltanRibBal = false;
  } else if( method == "graph" ) {
    parmaBal = false;
    zoltanBal = true;
    zoltanRibBal = false;
  } else if( method == "zrib" ) {
    parmaBal = false;
    zoltanBal = false;
    zoltanRibBal = true;
  } else if ( method == "none" ) {
    parmaBal = false;
    zoltanBal = false;
    zoltanRibBal = false;
  } else {
    if (!PCU_Comm_Self())
      fprintf(stderr,
          "warning: ignoring unknown value of %s = %s\n",
          key, method.c_str());
  }
}

void setupMatching(ma::Input& in) {
  if (!PCU_Comm_Self())
    printf("Matched mesh: disabling"
           " snapping, and shape correction,\n");
  in.shouldSnap = false;
  in.shouldFixShape = false;
}

struct AdaptCallback : public Parma_GroupCode
{
  apf::Mesh2* mesh;
  apf::Field* field;
  ph::Input* in;
  AdaptCallback(apf::Mesh2* m, apf::Field* szfld)
    : mesh(m), field(szfld), in(NULL) { }
  AdaptCallback(apf::Mesh2* m, apf::Field* szfld, ph::Input* inp)
    : mesh(m), field(szfld), in(inp) { }
  void run(int) {
    ma::Input* ma_in = ma::configure(mesh, field);
    if( in ) {
      //chef defaults
      ma_in->shouldRunPreZoltan = true;
      ma_in->shouldRunMidParma = true;
      ma_in->shouldRunPostParma = true;
      //override with user inputs if specified
      setupBalance("preAdaptBalanceMethod", in->preAdaptBalanceMethod,
          ma_in->shouldRunPreParma, ma_in->shouldRunPreZoltan,
          ma_in->shouldRunPreZoltanRib);
      bool ignored;
      setupBalance("midAdaptBalanceMethod", in->midAdaptBalanceMethod,
          ma_in->shouldRunMidParma, ma_in->shouldRunMidZoltan, ignored);
      setupBalance("postAdaptBalanceMethod", in->postAdaptBalanceMethod,
          ma_in->shouldRunPostParma, ma_in->shouldRunPostZoltan,
          ma_in->shouldRunPostZoltanRib);
      ma_in->shouldTransferParametric = in->transferParametric;
      ma_in->shouldSnap = in->snap;
      ma_in->maximumIterations = in->maxAdaptIterations;
      /*
        validQuality sets which elements will be accepted during mesh
        modification. If no boundary layers, you might bring this high (e.g,
        O(1e-2).  This would be bad for boundary layers though since their
        high aspect ratio routinely produces quality measures in the e-6 to e-7
        range so, when there are layers, this needs to be O(1e-8).
      */
      ma_in->validQuality = in->validQuality;
    } else {
      ma_in->shouldRunPreZoltan = true;
    }
    if (mesh->hasMatching())
      ph::setupMatching(*ma_in);
    ma::adapt(ma_in);
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
  PCU_ALWAYS_ASSERT(PCU_Comm_Peers() % factor == 0);
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
    fprintf(stderr,"adaptShrunken limit set to %f factor computed as %d\n", minPartDensity, factor);
  if (factor == 1) {
    callback.run(0);
  } else {
    warnAboutShrinking(factor);
    Parma_ShrinkPartition(m, factor, callback);
  }
}

static void runFromErrorThreshold(Input& in, apf::Mesh2* m)
{
  const char* fieldname = in.adaptErrorFieldName.c_str();
  const unsigned idx = in.adaptErrorFieldIndex;
  const double errLimit = in.adaptErrorThreshold;
  const double factor = 0.5;
  apf::Field* szFld = sam::errorThreshold(m,fieldname,idx,errLimit,factor);
  PCU_ALWAYS_ASSERT(szFld);
  chef::adapt(m, szFld);
  apf::destroyField(szFld);
}

static void runFromGivenSize(Input& in, apf::Mesh2* m)
{
  const unsigned idx = 5;
  apf::Field* szFld = sam::specifiedIso(m,"errors",idx);
  PCU_ALWAYS_ASSERT(szFld);
  chef::adapt(m, szFld,in);
  apf::destroyField(szFld);
}

void tetrahedronize(Input& in, apf::Mesh2* m)
{
  ma::Input* ma_in = ma::configureIdentity(m);
  ph::setupBalance("preAdaptBalanceMethod", in.preAdaptBalanceMethod,
      ma_in->shouldRunPreParma, ma_in->shouldRunPreZoltan,
      ma_in->shouldRunPreZoltanRib);
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
    ph::AdaptCallback acb(m,szFld);
    adaptShrunken(m,10000,acb);
  }

  void adapt(apf::Mesh2* m, apf::Field* szFld, ph::Input& in) {
    ph::AdaptCallback acb(m,szFld,&in);
    adaptShrunken(m,in.adaptShrinkLimit,acb);
  }

  void uniformRefinement(ph::Input& in, apf::Mesh2* m)
  {
    ma::Input* ma_in = ma::configureMatching(m, in.recursiveUR);
    ph::setupBalance("preAdaptBalanceMethod", in.preAdaptBalanceMethod,
        ma_in->shouldRunPreParma, ma_in->shouldRunPreZoltan,
        ma_in->shouldRunPreZoltanRib);
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

#include <stdio.h>
#include <PCU.h>
#include "parma.h"
#include "parma_balancer.h"
#include "parma_step.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"
#include "parma_commons.h"
#include "parma_convert.h"

namespace {
  using parmaCommons::status;

  class GhostElmBalancer : public parma::Balancer {
    private:
      int sideTol;
    public:
      GhostElmBalancer(apf::Mesh* m, int l, double f, int v)
        : Balancer(m, f, v, "ghostElms"), layers(l)
      {
        parma::Sides* s = parma::makeVtxSides(mesh);
        sideTol = TO_INT(parma::avgSharedSides(s));
        delete s;
        if( !PCU_Comm_Self() && verbose )
          status("sideTol %d\n", sideTol);
      }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        parma::Sides* s = parma::makeVtxSides(mesh);
        double avgSides = parma::avgSharedSides(s);
        if( !PCU_Comm_Self() && verbose )
          status("avgSides %f\n", avgSides);

        parma::GhostWeights* gw =
          parma::makeVtxGhostWeights(mesh, wtag, s, layers);
        parma::Weights* elmW = convertGhostToEntWeight(gw,mesh->getDimension());
        parma::Weights* vtxW =convertGhostToEntWeight(gw,0);
        destroyGhostWeights(gw);

        double vtxImb, vtxAvg;
        parma::getImbalance(vtxW, vtxImb, vtxAvg);
        if( !PCU_Comm_Self() && verbose )
          status("vtx imbalance %.3f avg %.3f\n", vtxImb, vtxAvg);
        delete vtxW;

        double elmImb, elmAvg;
        parma::getImbalance(elmW, elmImb, elmAvg);
        monitorUpdate(elmImb, iS, iA);
        monitorUpdate(avgSides, sS, sA);

        parma::Targets* t = parma::makeTargets(s, elmW, factor);
        parma::Selector* sel = parma::makeElmSelector(mesh, wtag);
        parma::BalOrStall* stopper =
          new parma::BalOrStall(iA, sA, sideTol*.001, verbose);
        parma::Stepper b(mesh, factor, s, elmW, t, sel, "elm", stopper);
        bool ret = b.step(tolerance, verbose);
        return ret;
      }
      int layers;
  };

  class GhostVtxLtElmBalancer : public parma::Balancer {
    private:
      int sideTol;
      int layers;
      double maxElmW;
      int stepNum;
    public:
      GhostVtxLtElmBalancer(apf::Mesh* m, double f, int v, int l)
        : Balancer(m, f, v, "ghostVtxLtElms"), layers(l), maxElmW(0), stepNum(0)
      {
        parma::Sides* s = parma::makeVtxSides(mesh);
        sideTol = TO_INT(parma::avgSharedSides(s));
        delete s;
        if( !PCU_Comm_Self() && verbose )
          status("sideTol %d\n", sideTol);
      }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        parma::Sides* s = parma::makeVtxSides(mesh);
        double avgSides = parma::avgSharedSides(s);
        if( !PCU_Comm_Self() && verbose )
          status("avgSides %f\n", avgSides);

        parma::GhostWeights* gw =
          parma::makeVtxGhostWeights(mesh, wtag, s, layers);
        parma::Weights* vtxW = convertGhostToEntWeight(gw,0);
        parma::Weights* edgeW = convertGhostToEntWeight(gw,1);
        parma::Weights* elmW = convertGhostToEntWeight(gw,mesh->getDimension());
        destroyGhostWeights(gw);

        double elmImb, elmAvg;
        parma::getImbalance(elmW,elmImb,elmAvg);
        double edgeImb, edgeAvg;
        parma::getImbalance(edgeW, edgeImb, edgeAvg);
        if( !PCU_Comm_Self() && verbose ) {
          status("elm imbalance %.3f avg %.3f\n", elmImb, elmAvg);
          status("edge imbalance %.3f avg %.3f\n", edgeImb, edgeAvg);
        }
        if( !stepNum ) //FIXME need to set the imbalance at the beginning for the primary entity
          maxElmW = parma::getMaxWeight(elmW);
        delete edgeW;

        monitorUpdate(elmImb, iS, iA);
        monitorUpdate(avgSides, sS, sA);
        parma::Targets* t =
          parma::makePreservingTargets(s, vtxW, elmW, sideTol, maxElmW, factor);
        delete elmW;
        parma::Selector* sel =
          parma::makeVtxLtElmSelector(mesh, wtag, maxElmW);
        parma::BalOrStall* stopper =
          new parma::BalOrStall(iA, sA, sideTol*.001, verbose);
        parma::Stepper b(mesh, factor, s, vtxW, t, sel, "vtx", stopper);
        bool ret = b.step(tolerance, verbose);
        stepNum++;
        return ret;
      }
  };
}

class GhostElmGtVtxBalancer : public parma::Balancer {
  public:
    GhostElmGtVtxBalancer(apf::Mesh* m, int l, double f, int v)
      : Balancer(m, f, v, "ghostsElmGtVtx"), layers(l) { }
    bool runStep(apf::MeshTag*, double) { return true; }
    void balance(apf::MeshTag* wtag, double tolerance) {
      apf::Balancer* b = new GhostElmBalancer(mesh, layers, factor, verbose);
      b->balance(wtag, tolerance);
      delete b;
      Parma_PrintPtnStats(mesh, "post-elements", (verbose>2));
      b = new GhostVtxLtElmBalancer(mesh, factor, verbose, layers);
      b->balance(wtag, tolerance);
      delete b;
    }
  private:
    int layers;
};

apf::Balancer* Parma_MakeGhostDiffuser(apf::Mesh* m,
    int layers, int, double stepFactor, int verbosity) {
  return Parma_MakeGhostDiffuser(m,layers,stepFactor,verbosity);
}

apf::Balancer* Parma_MakeGhostDiffuser(apf::Mesh* m,
    int layers, double stepFactor, int verbosity) {
  return new GhostElmGtVtxBalancer(m, layers, stepFactor, verbosity);
}

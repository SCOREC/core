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

  class GhostEdgeBalancer : public parma::Balancer {
    private:
      int sideTol;
    public:
      GhostEdgeBalancer(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "ghostEdges")
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
          parma::makeElmGhostWeights(mesh, wtag, s);
        parma::Weights* edgeW = convertGhostToEntWeight(gw,1);
        parma::Weights* faceW = convertGhostToEntWeight(gw,2);
        parma::Weights* elmW = convertGhostToEntWeight(gw,3);
        destroyGhostWeights(gw);

        double faceImb, faceAvg, elmImb, elmAvg;
        parma::getImbalance(faceW, faceImb, faceAvg);
        parma::getImbalance(elmW, elmImb, elmAvg);
        if( !PCU_Comm_Self() && verbose ) {
          status("face imbalance %.3f avg %.3f\n", faceImb, faceAvg);
          status("elm imbalance %.3f avg %.3f\n", elmImb, elmAvg);
        }
        delete faceW;
        delete elmW;

        double edgeImb, edgeAvg;
        parma::getImbalance(edgeW, edgeImb, edgeAvg);
        monitorUpdate(edgeImb, iS, iA);
        monitorUpdate(avgSides, sS, sA);

        parma::Targets* t = parma::makeTargets(s, edgeW, factor);
        parma::Selector* sel = parma::makeVtxSelector(mesh, wtag);
        parma::BalOrStall* stopper =
          new parma::BalOrStall(iA, sA, sideTol*.001, verbose);
        parma::Stepper b(mesh, factor, s, edgeW, t, sel, "edge", stopper);
        bool ret = b.step(tolerance, verbose);
        return ret;
      }
  };
}

apf::Balancer* Parma_MakeGhostEdgeDiffuser(apf::Mesh* m,
    double stepFactor, int verbosity) {
  return new GhostEdgeBalancer(m,stepFactor,verbosity);
}

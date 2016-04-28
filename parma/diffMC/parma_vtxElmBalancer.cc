#include <PCU.h>
#include <apfPartition.h>
#include <parma.h>
#include "parma_balancer.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"
#include "parma_step.h"
#include "parma_monitor.h"
#include "parma_commons.h"
#include "parma_convert.h"

namespace {
  using parmaCommons::status;

  class ElmLtVtx : public parma::Balancer {
    private:
      int sideTol;
      double maxVtx;
    public:
      ElmLtVtx(apf::Mesh* m, double f, double maxV, int v)
        : Balancer(m, f, v, "elements") {
          maxVtx = maxV;
          if( !PCU_Comm_Self() && verbose ) {
            status("stepFactor %.3f\n", f);
            status("maxVtx %.3f\n", maxVtx);
          }
          parma::Sides* s = parma::makeVtxSides(mesh);
          sideTol = TO_INT(parma::avgSharedSides(s));
          delete s;
          if( !PCU_Comm_Self() && verbose )
            status("sideTol %d\n", sideTol);
      }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        const double maxVtxImb =
          Parma_GetWeightedEntImbalance(mesh, wtag, 0);
        const double maxElmImb =
          Parma_GetWeightedEntImbalance(mesh, wtag, mesh->getDimension());
        if( !PCU_Comm_Self() && verbose )
          status("vtx imbalance %.3f\n", maxVtxImb);
        parma::Sides* s = parma::makeVtxSides(mesh);
        parma::Weights* vtxW = parma::makeEntWeights(mesh, wtag, s, 0);
        parma::Weights* elmW =
          parma::makeEntWeights(mesh, wtag, s, mesh->getDimension());
        parma::Targets* t =
          parma::makePreservingTargets(s, elmW, vtxW, sideTol, maxVtx, factor);
        delete vtxW;
        parma::Selector* sel =
          parma::makeElmLtVtxSelector(mesh, wtag, maxVtx);

        double avgSides = parma::avgSharedSides(s);
        monitorUpdate(maxElmImb, iS, iA);
        monitorUpdate(avgSides, sS, sA);
        if( !PCU_Comm_Self() && verbose )
          status("elmImb %f avgSides %f\n", maxElmImb, avgSides);
        parma::BalOrStall* stopper =
          new parma::BalOrStall(iA, sA, sideTol*.001, verbose);

        parma::Stepper b(mesh, factor, s, elmW, t, sel, "elm", stopper);
        return b.step(tolerance, verbose);
      }
  };
}

class VtxElmBalancer : public parma::Balancer {
  public:
    VtxElmBalancer(apf::Mesh* m, double f, int v)
      : Balancer(m, f, v, "cake") { }
    bool runStep(apf::MeshTag*, double) { return true; }
    void balance(apf::MeshTag* wtag, double tolerance) {
      apf::Balancer* b = Parma_MakeVtxBalancer(mesh, factor, verbose);
      b->balance(wtag, tolerance);
      delete b;
      Parma_PrintPtnStats(mesh, "post vertices", (verbose>2));
      double maxVtxW = parma::getMaxWeight(mesh, wtag, 0);
      b = new ElmLtVtx(mesh, factor, maxVtxW, verbose);
      b->balance(wtag, tolerance);
      delete b;
    }
};

apf::Balancer* Parma_MakeVtxElmBalancer(apf::Mesh* m,
    double stepFactor, int verbosity) {
  return new VtxElmBalancer(m, stepFactor, verbosity);
}

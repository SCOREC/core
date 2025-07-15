#include "parma.h"
#include "parma_step.h"
#include "parma_balancer.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"
#include "parma_monitor.h"
#include "parma_stop.h"
#include "parma_graphDist.h"
#include "parma_commons.h"
#include "parma_convert.h"

namespace {
  using parmaCommons::status;

  class VtxBalancer : public parma::Balancer {
    private:
      int sideTol;
    public:
      VtxBalancer(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "vertices") {
          parma::Sides* s = parma::makeVtxSides(mesh);
          sideTol = TO_INT(parma::avgSharedSides(s, mesh->getPCU()));
          delete s;
          if( !m->getPCU()->Self() && verbose )
            status("sideTol %d\n", sideTol);
      }

      bool runStep(apf::MeshTag* wtag, double tolerance) {
        const double maxVtxImb =
          Parma_GetWeightedEntImbalance(mesh, wtag, 0);
        parma::Sides* s = parma::makeVtxSides(mesh);
        parma::Weights* w = parma::makeEntWeights(mesh, wtag, s, 0);
        parma::Targets* t =
          parma::makeWeightSideTargets(s, w, sideTol, factor);
        parma::Selector* sel = parma::makeVtxSelector(mesh, wtag);
        double avgSides = parma::avgSharedSides(s, mesh->getPCU());
        monitorUpdate(maxVtxImb, iS, iA);
        monitorUpdate(avgSides, sS, sA);
        if( !mesh->getPCU()->Self() && verbose )
          status("vtxImb %f avgSides %f\n", maxVtxImb, avgSides);
        parma::BalOrStall* stopper = 
          new parma::BalOrStall(iA, sA, sideTol*.001, verbose);
        parma::Stepper b(mesh, factor, s, w, t, sel, "vtx", stopper);
        return b.step(tolerance, verbose);
      }
  };
}

apf::Balancer* Parma_MakeVtxBalancer(apf::Mesh* m,
    double stepFactor, int verbosity) {
  if( !m->getPCU()->Self() && verbosity )
    status("stepFactor %.3f\n", stepFactor);
  return new VtxBalancer(m, stepFactor, verbosity);
}

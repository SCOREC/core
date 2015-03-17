#include <PCU.h>
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

namespace {
  class VtxBalancer : public parma::Balancer {
    private:
      int sideTol;
    public:
      VtxBalancer(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "vertices") {
          parma::Sides* s = parma::makeVtxSides(mesh);
          sideTol = static_cast<int>(parma::avgSharedSides(s));
          delete s;
          if( !PCU_Comm_Self() )
            fprintf(stdout, "sideTol %d\n", sideTol);
      }

      bool runStep(apf::MeshTag* wtag, double tolerance) {
        const double maxVtxImb =
          Parma_GetWeightedEntImbalance(mesh, wtag, 0);
        parma::Sides* s = parma::makeVtxSides(mesh);
        parma::Weights* w = parma::makeEntWeights(mesh, wtag, s, 0);
        parma::Targets* t =
          parma::makeWeightSideTargets(s, w, sideTol, factor);
        parma::Selector* sel = parma::makeVtxSelector(mesh, wtag);
        double avgSides = parma::avgSharedSides(s);
        monitorUpdate(maxVtxImb, iS, iA);
        monitorUpdate(avgSides, sS, sA);
        if( !PCU_Comm_Self() )
          fprintf(stdout, "vtxImb %f avgSides %f\n", maxVtxImb, avgSides);
        parma::BalOrStall* stopper = 
          new parma::BalOrStall(iA, sA, sideTol*.001);
        parma::Stepper b(mesh, factor, s, w, t, sel, stopper);
        return b.step(tolerance, verbose);
      }
  };
}

apf::Balancer* Parma_MakeVtxBalancer(apf::Mesh* m,
    double stepFactor, int verbosity) {
  if( !PCU_Comm_Self() && verbosity )
    fprintf(stdout,"PARMA_STATUS stepFactor %.3f\n", stepFactor);
  return new VtxBalancer(m, stepFactor, verbosity);
}

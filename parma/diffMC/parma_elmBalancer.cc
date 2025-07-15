#include <parma_balancer.h>
#include "parma.h"
#include "parma_step.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"
#include "parma_commons.h"

namespace {
  using parmaCommons::status;

  class ElmBalancer : public parma::Balancer {
    private:
      double sideTol;
    public:
      ElmBalancer(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "elements") {
          parma::Sides* s = parma::makeVtxSides(mesh);
          sideTol = parma::avgSharedSides(s, mesh->getPCU());
          delete s;
      }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        const double maxElmImb =
          Parma_GetWeightedEntImbalance(mesh, wtag, mesh->getDimension());
        parma::Sides* s = parma::makeVtxSides(mesh);
        double avgSides = parma::avgSharedSides(s, mesh->getPCU());
        parma::Weights* w =
          parma::makeEntWeights(mesh, wtag, s, mesh->getDimension());
        parma::Targets* t = parma::makeTargets(s, w, factor);
        parma::Selector* sel = parma::makeElmSelector(mesh, wtag);

        monitorUpdate(maxElmImb, iS, iA);
        monitorUpdate(avgSides, sS, sA);
        if( !mesh->getPCU()->Self() && verbose )
          status("elmImb %f avgSides %f\n", maxElmImb, avgSides);
        parma::BalOrStall* stopper =
          new parma::BalOrStall(iA, sA, sideTol*.001, verbose);

        parma::Stepper b(mesh, factor, s, w, t, sel, "elm", stopper);
        return b.step(tolerance, verbose);
      }
  };
}

apf::Balancer* Parma_MakeElmBalancer(apf::Mesh* m,
    double stepFactor, int verbosity) {
  if( !m->getPCU()->Self() && verbosity )
    status("stepFactor %.3f\n", stepFactor);
  return new ElmBalancer(m, stepFactor, verbosity);
}

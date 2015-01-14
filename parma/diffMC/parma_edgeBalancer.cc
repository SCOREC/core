#include "parma.h"
#include "parma_balancer.h"
#include "parma_step.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"

namespace {
  class EdgeBalancer : public parma::Balancer {
    public:
      EdgeBalancer(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "edges") { }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        parma::Sides* s = parma::makeVtxSides(mesh);
        parma::Weights* w = parma::makeEntWeights(mesh, wtag, s, 1);
        parma::Targets* t = parma::makeTargets(s, w, factor);
        parma::Selector* sel = parma::makeEdgeSelector(mesh, wtag);
        parma::Stepper b(mesh, factor, s, w, t, sel);
        return b.step(tolerance, verbose);
      }
  };
}

apf::Balancer* Parma_MakeEdgeBalancer(apf::Mesh* m, 
    double stepFactor, int verbosity) {
  return new EdgeBalancer(m, stepFactor, verbosity);
}

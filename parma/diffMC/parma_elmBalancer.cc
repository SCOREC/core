#include <PCU.h>
#include <parma_balancer.h>
#include "parma_step.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"

namespace {
  class ElmBalancer : public parma::Balancer {
    public:
      ElmBalancer(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "elements") { }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        parma::Sides* s = parma::makeElmBdrySides(mesh);
        parma::Weights* w =
          parma::makeEntWeights(mesh, wtag, s, mesh->getDimension());
        parma::Targets* t = parma::makeTargets(s, w, factor);
        parma::Selector* sel = parma::makeElmSelector(mesh, wtag);
        parma::Stepper b(mesh, wtag, factor, s, w, t, sel);
        return b.step(tolerance, verbose);
      }
  };
}

apf::Balancer* Parma_MakeElmBalancer(apf::Mesh* m,
    double stepFactor, int verbosity) {
  return new ElmBalancer(m, stepFactor, verbosity);
}

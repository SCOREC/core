#include <stdio.h>
#include "parma.h"
#include "parma_balancer.h"
#include "parma_step.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"

namespace {
  class GhostBalancer : public parma::Balancer {
    public:
      GhostBalancer(apf::Mesh* m, int l, int b, double f, int v)
        : Balancer(m, f, v, "ghosts"), layers(l), bridge(b) { }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        parma::Sides* s = parma::makeElmBdrySides(mesh);
        parma::Weights* w =
          parma::makeGhostWeights(mesh, wtag, s, layers, bridge);
        parma::Targets* t = parma::makeTargets(s, w, factor);
        parma::Selector* sel = parma::makeVtxSelector(mesh, wtag);
        parma::Stepper b(mesh, factor, s, w, t, sel);
        bool ret = b.step(tolerance, verbose);
        return ret;
      }
      int layers;
      int bridge;
  };
}

apf::Balancer* Parma_MakeGhostDiffuser(apf::Mesh* m,
    int layers, int bridge, double stepFactor, int verbosity) {
  return new GhostBalancer(m, layers, bridge, stepFactor, verbosity);
}

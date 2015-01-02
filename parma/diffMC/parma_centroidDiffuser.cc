#include <PCU.h>
#include <stdio.h>
#include "parma.h"
#include "parma_step.h"
#include "parma_balancer.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_centroids.h"
#include "parma_targets.h"
#include "parma_selector.h"

namespace {
  class CentroidBalancer : public parma::Balancer {
    public:
      CentroidBalancer(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "elements") { }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        parma::Sides* s = parma::makeElmBdrySides(mesh);
        parma::Weights* w =
          parma::makeEntWeights(mesh, wtag, s, mesh->getDimension());
        parma::Targets* t = parma::makeTargets(s, w, factor);
        parma::Centroids c(mesh, wtag, s);
        parma::Selector* sel = parma::makeCentroidSelector(mesh, wtag, &c);
        parma::Stepper b(mesh, factor, s, w, t, sel);
        return b.step(tolerance, verbose);
      }
  };
}

apf::Balancer* Parma_MakeCentroidDiffuser(apf::Mesh* m,
    double stepFactor, int verbosity) {
  return new CentroidBalancer(m, stepFactor, verbosity);
}

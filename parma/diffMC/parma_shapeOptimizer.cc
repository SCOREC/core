#include <PCU.h>
#include <stdio.h>
#include "parma_balancer.h"
#include "parma_step.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_centroids.h"
#include "parma_targets.h"
#include "parma_selector.h"

namespace {
  class ShapeOptimizer : public parma::Balancer {
    public:
      ShapeOptimizer(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "gap") { }
      static bool greater(double imb, double maxImb) {
        return imb > maxImb;
      }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        parma::Sides* s = parma::makeElmBdrySides(mesh);
        parma::Weights* w =
          parma::makeEntWeights(mesh, wtag, s, mesh->getDimension());
        parma::Targets* t = parma::makeShapeTargets(mesh, s, w, factor);
        PCU_Debug_Print("%s\n", t->print("targets").c_str());
        parma::Centroids c(mesh, wtag, s);
        parma::Selector* sel = parma::makeShapeSelector(mesh, wtag, &c);
        parma::Stepper b(mesh, wtag, factor, s, w, t, sel, greater);
        return b.step(tolerance, verbose);
      }
  };
}

apf::Balancer* Parma_MakeShapeOptimizer(apf::Mesh* m,
    double stepFactor, int verbosity) {
  return new ShapeOptimizer(m, stepFactor, verbosity);
}

#include <apfPartition.h>
#include <parma.h>
#include "parma_balancer.h"
#include "parma_step.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"


namespace {
  class VtxEdgeBalancer : public parma::Balancer {
    public:
      VtxEdgeBalancer(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "edges") { }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        parma::Sides* s = parma::makeElmBdrySides(mesh);
        parma::Weights* w[2] =
          {parma::makeEntWeights(mesh, wtag, s, 0),
            parma::makeEntWeights(mesh, wtag, s, 1)};
        parma::Targets* t = parma::makeVtxEdgeTargets(s, w, factor);
        parma::Selector* sel = parma::makeEdgeSelector(mesh, wtag);
        parma::Stepper b(mesh, factor, s, w[1], t, sel);
        return b.step(tolerance, verbose);
      }
  };

  class VtxEdgeElmBalancer : public parma::Balancer {
    public:
      VtxEdgeElmBalancer(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "cake") { }
      bool runStep(apf::MeshTag*, double) { return true; }
      void balance(apf::MeshTag* wtag, double tolerance) {
        apf::Balancer* b = Parma_MakeElmBalancer(mesh, factor, verbose);
        b->balance(wtag, tolerance);
        Parma_PrintPtnStats(mesh, "post elements");
        delete b;

        b = Parma_MakeVtxBalancer(mesh, factor, verbose);
        b->balance(wtag, tolerance);
        Parma_PrintPtnStats(mesh, "post vertices");
        delete b;

        b = new VtxEdgeBalancer(mesh, factor, verbose);
        b->balance(wtag, tolerance);
        Parma_PrintPtnStats(mesh, "post edges");
        delete b;
      }
  };
}

apf::Balancer* Parma_MakeVtxEdgeElmBalancer(apf::Mesh* m,
    double stepFactor, int verbosity) {
  return new VtxEdgeElmBalancer(m, stepFactor, verbosity);
}

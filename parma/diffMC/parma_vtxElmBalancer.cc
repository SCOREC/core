#include <apfPartition.h>
#include <parma.h>
#include "parma_balancer.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"
#include "parma_step.h"

namespace {
  class ElmLtVtx : public parma::Balancer {
    public:
      ElmLtVtx(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "elements") { }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        parma::Sides* s = parma::makeElmBdrySides(mesh);
        parma::Weights* w[2] =
        {parma::makeEntWeights(mesh, wtag, s, 0),
          parma::makeEntWeights(mesh, wtag, s, mesh->getDimension())};
        parma::Targets* t = parma::makeVtxElmTargets(s, w, factor);
        parma::Selector* sel = parma::makeElmSelector(mesh, wtag);
        parma::Stepper b(mesh, wtag, factor, s, w[1], t, sel);
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
        Parma_PrintPtnStats(mesh, "post vertices");
        b = new ElmLtVtx(mesh, factor, verbose);
        b->balance(wtag, tolerance);
        delete b;
      }
  };


apf::Balancer* Parma_MakeVtxElmBalancer(apf::Mesh* m,
    double stepFactor, int verbosity) {
  return new VtxElmBalancer(m, stepFactor, verbosity);
}

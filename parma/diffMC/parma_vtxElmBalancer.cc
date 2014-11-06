#include <apfPartition.h>
#include <parma.h>
#include "parma_balancer.h"

namespace {
  class VtxElmBalancer : public parma::Balancer {
    public:
      VtxElmBalancer(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "cake") { }
      bool runStep(apf::MeshTag*, double) { return true; }
      void balance(apf::MeshTag* wtag, double tolerance) {
        apf::Balancer* b = Parma_MakeElmBalancer(mesh, factor, verbose);
        b->balance(wtag, tolerance);
        delete b;
        b = Parma_MakeVtxBalancer(mesh, factor, verbose);
        b->balance(wtag, tolerance);
        delete b;
      }
  };
}

apf::Balancer* Parma_MakeVtxElmBalancer(apf::Mesh* m,
    double stepFactor, int verbosity) {
  return new VtxElmBalancer(m, stepFactor, verbosity);
}

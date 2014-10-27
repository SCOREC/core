#include <PCU.h>
#include <apfPartition.h>
#include "parma_base.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"

namespace {
  void printTiming(const char* type, int steps, double tol, double time) {
    if (!PCU_Comm_Self())
      printf("%s balanced in %d steps to %f in %f seconds\n",
          type, steps, tol, time);
  }
}

namespace parma {
  class ElmBalancer : public apf::Balancer {
    public:
      ElmBalancer(apf::Mesh* m, double f, int v)
        : mesh(m), factor(f), verbose(v) {
          (void) verbose; // silence!
        }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        Sides* s = makeElmBdrySides(mesh);
        Weights* w = makeEntWeights(mesh, wtag, s, mesh->getDimension());
        Targets* t = makeTargets(s, w, factor);
        Selector* sel = makeElmSelector(mesh, wtag);
        parma::Balancer b(mesh, wtag, factor, s, w, t, sel);
        return b.run(tolerance, verbose);
      }
      virtual void balance(apf::MeshTag* wtag, double tolerance) {
        int step = 0;
        double t0 = MPI_Wtime();
        while (runStep(wtag,tolerance) && step++ < 50);
        printTiming("elements", step, tolerance, MPI_Wtime()-t0);
      }
    private:
      apf::Mesh* mesh;
      double factor;
      int verbose;
  };
} //end parma namespace

apf::Balancer* Parma_MakeElmBalancer(apf::Mesh* m, 
    double stepFactor, int verbosity) {
  return new parma::ElmBalancer(m, stepFactor, verbosity);
}

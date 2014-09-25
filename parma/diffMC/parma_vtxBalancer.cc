#include <PCU.h>
#include <apfPartition.h>
#include "parma_base.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"

namespace parma {
  class VtxBalancer : public apf::Balancer {
    public:
      VtxBalancer(apf::Mesh* m, double f, int v)
        : mesh(m), factor(f), verbose(v) {
          (void) verbose; // silence!
        }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        Sides* s = makeElmBdrySides(mesh);
        Weights* w = makeEntWeights(mesh, wtag, s, 0);
        Targets* t = makeTargets(s, w, factor);
        Selector* sel = makeVtxSelector(mesh, wtag);
        parma::Balancer parmaVtx(mesh, wtag, factor, s, w, t, sel);
        return parmaVtx.run(tolerance);
      }
      virtual void balance(apf::MeshTag* wtag, double tolerance) {
        double t0 = MPI_Wtime();
        while (runStep(wtag,tolerance));
        double t1 = MPI_Wtime();
        if (!PCU_Comm_Self())
          printf("vertices balanced to %f in %f seconds\n", tolerance, t1-t0);
      }
    private:
      apf::Mesh* mesh;
      double factor;
      int verbose;
  };
} //end parma namespace

apf::Balancer* Parma_MakeVtxBalancer(apf::Mesh* m, 
    double stepFactor, int verbosity) {
  return new parma::VtxBalancer(m, stepFactor, verbosity);
}

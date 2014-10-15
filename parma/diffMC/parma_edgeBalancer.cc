#include <PCU.h>
#include <apfPartition.h>
#include "parma_base.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"

namespace parma {
  class EdgeBalancer : public apf::Balancer {
    public:
      EdgeBalancer(apf::Mesh* m, double f, int v)
        : mesh(m), factor(f), verbose(v) {
          (void) verbose; // silence!
        }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        Sides* s = makeElmBdrySides(mesh);
        Weights* w = makeEntWeights(mesh, wtag, s, 1);
        Targets* t = makeTargets(s, w, factor);
        Selector* sel = makeEdgeSelector(mesh, wtag);
        parma::Balancer parmaEdge(mesh, wtag, factor, s, w, t, sel);
        return parmaEdge.run(tolerance, 1);
      }
      virtual void balance(apf::MeshTag* wtag, double tolerance) {
        double t0 = MPI_Wtime();
        int steps = 0;
        while (runStep(wtag,tolerance) && steps++ < 100);
        double t1 = MPI_Wtime();
        if (!PCU_Comm_Self())
          printf("edges balanced to %f in %f seconds\n", tolerance, t1-t0);
      }
    private:
      apf::Mesh* mesh;
      double factor;
      int verbose;
  };
} //end parma namespace

apf::Balancer* Parma_MakeEdgeBalancer(apf::Mesh* m, 
    double stepFactor, int verbosity) {
  return new parma::EdgeBalancer(m, stepFactor, verbosity);
}

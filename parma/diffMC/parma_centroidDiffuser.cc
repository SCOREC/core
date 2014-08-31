#include <PCU.h>
#include <stdio.h>
#include "parma.h"
#include "parma_base.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_centroids.h"
#include "parma_targets.h"
#include "parma_selector.h"

namespace parma {
  class CentroidBalancer : public apf::Balancer {
    public:
      CentroidBalancer(apf::Mesh* m, double f, int v)
        : mesh(m), factor(f), verbose(v) {
          (void) verbose; // silence!
      }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        Sides* s = makeElmBdrySides(mesh);
        Weights* w = makeEntWeights(mesh, wtag, s, mesh->getDimension());
        Targets* t = makeTargets(s, w, factor);
        Centroids c(mesh, wtag, s);
        Selector* sel = makeCentroidSelector(mesh, wtag, &c);
        parma::Balancer b(mesh, wtag, factor, s, w, t, sel); 
        return b.run(tolerance);
      }
      void balance(apf::MeshTag* weights, double tolerance) {
        double t0 = MPI_Wtime(); 
        while (runStep(weights,tolerance));
        double t1 = MPI_Wtime();
        if (!PCU_Comm_Self())
          printf("centroid balanced to %f in %f seconds\n", tolerance, t1-t0);
      }
    private:
      apf::Mesh* mesh;
      double factor;
      int verbose;
  };
}

apf::Balancer* Parma_MakeCentroidDiffuser(apf::Mesh* m, 
    double stepFactor, int verbosity) {
  return new parma::CentroidBalancer(m, stepFactor, verbosity);
}

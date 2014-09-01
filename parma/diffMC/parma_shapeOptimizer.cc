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
  class ShapeOptimizer : public apf::Balancer {
    public:
      ShapeOptimizer(apf::Mesh* m, double f, int v)
        : mesh(m), factor(f), verbose(v) {
          (void) verbose; // silence!
      }
      static bool greater(double imb, double maxImb) {
        return imb > maxImb;
      }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        Sides* s = makeElmBdrySides(mesh);
        Weights* w = makeEntWeights(mesh, wtag, s, mesh->getDimension());
        Targets* t = makeShapeTargets(mesh, s, w, factor);
        PCU_Debug_Print("%s\n", t->print("targets").c_str());
        Centroids c(mesh, wtag, s);
        Selector* sel = makeShapeSelector(mesh, wtag, &c);
        parma::Balancer b(mesh, wtag, factor, s, w, t, sel, greater); 
        return b.run(tolerance, verbose);
      }
      void balance(apf::MeshTag* weights, double tolerance) {
        double t0 = MPI_Wtime(); 
        int step = 0;
        while (runStep(weights,tolerance) && ++step );
        double t1 = MPI_Wtime();
        if (!PCU_Comm_Self())
          printf("gap ran %d steps in %f seconds\n", step, t1-t0);
      }
    private:
      apf::Mesh* mesh;
      double factor;
      int verbose;
  };
}

apf::Balancer* Parma_MakeShapeOptimizer(apf::Mesh* m, 
    double stepFactor, int verbosity) {
  return new parma::ShapeOptimizer(m, stepFactor, verbosity);
}

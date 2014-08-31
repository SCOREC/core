#include <PCU.h>
#include <stdio.h>
#include "parma.h"
#include "parma_base.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"

namespace parma {
  class GhostBalancer : public apf::Balancer {
    public:
      GhostBalancer(apf::Mesh* m, int l, int b, double f, int v)
        : mesh(m), factor(f), layers(l), bridge(b), verbose(v) {
          (void) verbose; // silence!
        }
      void printStats(Weights* w) {
        double max = w->self();
        PCU_Max_Doubles(&max, 1);
        double avg = w->self();
        PCU_Add_Doubles(&avg, 1);
        if ( !PCU_Comm_Self() )
          fprintf(stdout, "Max weight %.3f Avg weight %.3f\n", max, avg/PCU_Comm_Peers());
      }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        Sides* s = makeElmBdrySides(mesh);
        Weights* w = makeGhostWeights(mesh, wtag, s, layers, bridge);
        printStats(w);
        Targets* t = makeTargets(s, w, factor);
        Selector* sel = makeVtxSelector(mesh, wtag);
        parma::Balancer ghost(mesh, wtag, factor, s, w, t, sel);
        bool ret = ghost.run(tolerance);
        return ret;
      }
      void balance(apf::MeshTag* weights, double tolerance) {
        double t0 = MPI_Wtime(); 
        while (runStep(weights,tolerance));
        double t1 = MPI_Wtime();
        if (!PCU_Comm_Self())
          printf("ghost balanced to %f in %f seconds\n", tolerance, t1-t0);
      }
    private:
      apf::Mesh* mesh;
      double factor;
      int layers;
      int bridge;
      int verbose;
  };
}

apf::Balancer* Parma_MakeGhostDiffuser(apf::Mesh* m, 
    int layers, int bridge, double stepFactor, int verbosity) {
  return new parma::GhostBalancer(m, layers, bridge, stepFactor, verbosity);
}

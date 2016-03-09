#include <stdio.h>
#include <PCU.h>
#include "parma.h"
#include "parma_balancer.h"
#include "parma_step.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"

namespace {
  class GhostElementBalancer : public parma::Balancer {
    private: 
      int sideTol;
    public:
      GhostElementBalancer(apf::Mesh* m, int l, int b, double f, int v)
        : Balancer(m, f, v, "ghosts"), layers(l), bridge(b) 
      {
        parma::Sides* s = parma::makeElmBdrySides(mesh);
        sideTol = static_cast<int>(parma::avgSharedSides(s));
        delete s;
        if( !PCU_Comm_Self() && verbose )
          fprintf(stdout, "sideTol %d\n", sideTol);
      }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        parma::Sides* s = parma::makeElmBdrySides(mesh);

        const double maxVtxImb =
          Parma_GetWeightedEntImbalance(mesh, wtag, 0);
        double avgSides = parma::avgSharedSides(s);
        monitorUpdate(maxVtxImb, iS, iA);
        monitorUpdate(avgSides, sS, sA);
        if( !PCU_Comm_Self() && verbose )
          fprintf(stdout, "avgSides %f\n", avgSides);

        parma::Weights* w =
          parma::makeGhostElementWeights(mesh, wtag, s, layers, bridge);
        parma::Targets* t = parma::makeTargets(s, w, factor);
        parma::Selector* sel = parma::makeVtxSelector(mesh, wtag);

        /*parma::BalOrStall* stopper = 
          new parma::BalOrStall(iA, sA, sideTol*.001, verbose);
	*/
	parma::Stepper b(mesh, factor, s, w, t, sel, "elm");
        bool ret = b.step(tolerance, verbose);
        return ret;
      }
      int layers;
      int bridge;
  };
}

apf::Balancer* Parma_MakeGhostElementDiffuser(apf::Mesh* m,
    int layers, int bridge, double stepFactor, int verbosity) {
  return new GhostElementBalancer(m, layers, bridge, stepFactor, verbosity);
}

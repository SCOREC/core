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
  class MPASGhostBalancer : public parma::Balancer {
    private: 
      int sideTol;
    public:
      MPASGhostBalancer(apf::Mesh* m, int l, int b, double f, int v)
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

        const double maxElmImb =
          Parma_GetWeightedEntImbalance(mesh, wtag, mesh->getDimension());
        double avgSides = parma::avgSharedSides(s);
        monitorUpdate(maxElmImb, iS, iA);
        monitorUpdate(avgSides, sS, sA);
        if( !PCU_Comm_Self() && verbose )
          fprintf(stdout, "avgSides %f\n", avgSides);

        parma::Weights* w =
          parma::makeGhostMPASWeights(mesh, wtag, s, layers, bridge);
        parma::Targets* t = parma::makeTargets(s, w, factor);
        parma::Selector* sel = parma::makeVtxSelector(mesh, wtag);

        parma::BalOrStall* stopper =
          new parma::BalOrStall(iA, sA, sideTol*.001, verbose);
	parma::Stepper b(mesh, factor, s, w, t, sel, "elm", stopper);
        bool ret = b.step(tolerance, verbose);
        return ret;
      }
      int layers;
      int bridge;
  };
}

apf::Balancer* Parma_MakeMPASDiffuser(apf::Mesh* m,
    int layers, int bridge, double stepFactor, int verbosity) {
  return new MPASGhostBalancer(m, layers, bridge, stepFactor, verbosity);
}

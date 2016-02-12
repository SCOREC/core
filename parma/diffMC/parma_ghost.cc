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
  class GhostElmBalancer : public parma::Balancer {
    private: 
      int sideTol;
    public:
      GhostElmBalancer(apf::Mesh* m, int l, double f, int v)
        : Balancer(m, f, v, "ghostElms"), layers(l)
      {
        parma::Sides* s = parma::makeVtxSides(mesh);
        sideTol = static_cast<int>(parma::avgSharedSides(s));
        delete s;
        if( !PCU_Comm_Self() && verbose )
          fprintf(stdout, "sideTol %d\n", sideTol);
      }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        parma::Sides* s = parma::makeVtxSides(mesh);

        const double maxElmImb =
          Parma_GetWeightedEntImbalance(mesh, wtag, mesh->getDimension());
        double avgSides = parma::avgSharedSides(s);
        monitorUpdate(maxElmImb, iS, iA);
        monitorUpdate(avgSides, sS, sA);
        if( !PCU_Comm_Self() && verbose )
          fprintf(stdout, "avgSides %f\n", avgSides);

        parma::GhostWeights* gw =
          parma::makeGhostWeights(mesh, wtag, s, layers);
        parma::Weights* elmW = convertGhostToEntWeight(gw,3);
        destroyGhostWeights(gw);

        parma::Targets* t = parma::makeTargets(s, elmW, factor);
        parma::Selector* sel = parma::makeVtxSelector(mesh, wtag);

        if( verbose > 3 ) {
          fprintf(stderr, "%d tot %d %s\n", PCU_Comm_Self(), s->total(), s->print("sides").c_str());
          fprintf(stderr, "%d self %f %s\n", PCU_Comm_Self(), elmW->self(), elmW->print("weights").c_str());
          fprintf(stderr, "%d %s\n", PCU_Comm_Self(), t->print("tgts").c_str());
        }

        parma::BalOrStall* stopper = 
          new parma::BalOrStall(iA, sA, sideTol*.001, verbose);
        parma::Stepper b(mesh, factor, s, elmW, t, sel, stopper);
        bool ret = b.step(tolerance, verbose);
        return ret;
      }
      int layers;
  };
}

class GhostElmGtVtxBalancer : public parma::Balancer {
  public:
    GhostElmGtVtxBalancer(apf::Mesh* m, int l, double f, int v)
      : Balancer(m, f, v, "ghostsElmGtVtx"), layers(l) { }
    bool runStep(apf::MeshTag*, double) { return true; }
    void balance(apf::MeshTag* wtag, double tolerance) {
      apf::Balancer* b = new GhostElmBalancer(mesh, layers, factor, verbose);
      b->balance(wtag, tolerance);
      delete b;
    }
  private:
    int layers;
};


apf::Balancer* Parma_MakeGhostDiffuser(apf::Mesh* m,
    int layers, double stepFactor, int verbosity) {
  return new GhostElmGtVtxBalancer(m, layers, stepFactor, verbosity);
}

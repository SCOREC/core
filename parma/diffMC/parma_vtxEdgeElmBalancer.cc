#include <PCU.h>
#include <apfPartition.h>
#include <parma.h>
#include "parma_balancer.h"
#include "parma_step.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"


namespace {
  class VtxEdgeBalancer : public parma::Balancer {
    private:
      int sideTol;
      double maxVtx;
    public:
      VtxEdgeBalancer(apf::Mesh* m, double f, double maxV, int v)
        : Balancer(m, f, v, "edges") {
          maxVtx = maxV;
          if( !PCU_Comm_Self() ) {
            fprintf(stdout, "PARMA_STATUS stepFactor %.3f\n", f);
            fprintf(stdout, "PARMA_STATUS maxVtx %.3f\n", maxVtx);
          }
          parma::Sides* s = parma::makeVtxSides(mesh);
          sideTol = static_cast<int>(parma::avgSharedSides(s));
          delete s;
          if( !PCU_Comm_Self() )
            fprintf(stdout, "sideTol %d\n", sideTol);
      }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        const double maxVtxImb =
          Parma_GetWeightedEntImbalance(mesh, wtag, 0);
        if( !PCU_Comm_Self() )
          fprintf(stdout, "vtx imbalance %.3f\n", maxVtxImb);
        const double maxEdgeImb =
          Parma_GetWeightedEntImbalance(mesh, wtag, 1);
        parma::Sides* s = parma::makeVtxSides(mesh);
        double avgSides = parma::avgSharedSides(s);
        parma::Weights* w[2] =
          {parma::makeEntWeights(mesh, wtag, s, 0),
            parma::makeEntWeights(mesh, wtag, s, 1)};
        parma::Targets* t =
          parma::makeVtxEdgeTargets(s, w, sideTol, maxVtx, factor);
        parma::Selector* sel = parma::makeEdgeSelector(mesh, wtag);

        monitorUpdate(maxEdgeImb, iS, iA);
        monitorUpdate(avgSides, sS, sA);
        if( !PCU_Comm_Self() )
          fprintf(stdout, "edgeImb %f avgSides %f\n", maxEdgeImb, avgSides);
        parma::BalOrStall* stopper =
          new parma::BalOrStall(iA, sA, sideTol*.001);

        parma::Stepper b(mesh, factor, s, w[1], t, sel, stopper);
        bool ok = b.step(tolerance, verbose);
        delete w[0];
        return ok;
      }
  };

  class VtxEdgeElmBalancer : public parma::Balancer {
    public:
      VtxEdgeElmBalancer(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "cake") { }
      bool runStep(apf::MeshTag*, double) { return true; }
      void balance(apf::MeshTag* wtag, double tolerance) {
        apf::Balancer* b = Parma_MakeVtxBalancer(mesh, factor, verbose);
        b->balance(wtag, tolerance);
        Parma_PrintPtnStats(mesh, "post vertices");
        delete b;

        double maxVtxW = parma::getMaxWeight(mesh, wtag, 0);
        b = new VtxEdgeBalancer(mesh, factor, maxVtxW, verbose);
        b->balance(wtag, tolerance);
        Parma_PrintPtnStats(mesh, "post edges");
        delete b;

        maxVtxW = parma::getMaxWeight(mesh, wtag, 0);
        b = Parma_MakeElmLtVtxBalancer(mesh, maxVtxW, factor, verbose);
        b->balance(wtag, tolerance);
        Parma_PrintPtnStats(mesh, "post elements");
        delete b;
      }
  };
}

apf::Balancer* Parma_MakeVtxEdgeElmBalancer(apf::Mesh* m,
    double stepFactor, int verbosity) {
  return new VtxEdgeElmBalancer(m, stepFactor, verbosity);
}

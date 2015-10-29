#include <PCU.h>
#include <apfPartition.h>
#include <parma.h>
#include "parma_balancer.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"
#include "parma_step.h"
#include "parma_monitor.h"

namespace {
  class ElmLtVtxEdge : public parma::Balancer {
    private:
      int sideTol;
      double maxVtx;
      double maxEdge;
    public:
      ElmLtVtxEdge(apf::Mesh* m, double f, double maxV, double maxE, int v)
        : Balancer(m, f, v, "elements") {
          maxVtx = maxV;
          maxEdge = maxE;
          if( !PCU_Comm_Self() && verbose ) {
            fprintf(stdout, "PARMA_STATUS stepFactor %.3f\n", f);
            fprintf(stdout, "PARMA_STATUS maxVtx %.3f\n", maxVtx);
            fprintf(stdout, "PARMA_STATUS maxEdge %.3f\n", maxEdge);
          }
          parma::Sides* s = parma::makeVtxSides(mesh);
          sideTol = static_cast<int>(parma::avgSharedSides(s));
          delete s;
          if( !PCU_Comm_Self() && verbose )
            fprintf(stdout, "sideTol %d\n", sideTol);
      }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        const double maxElmImb =
          Parma_GetWeightedEntImbalance(mesh, wtag, mesh->getDimension());
        parma::Sides* s = parma::makeVtxSides(mesh);
        parma::Weights* w[3] =
          {parma::makeEntWeights(mesh, wtag, s, 0),
           parma::makeEntWeights(mesh, wtag, s, 1),
           parma::makeEntWeights(mesh, wtag, s, mesh->getDimension())};
        parma::Targets* t =
          parma::makeElmLtVtxEdgeTargets(s, w, sideTol, maxVtx, maxEdge, factor);
        delete w[0];
        delete w[1];
        parma::Selector* sel =
          parma::makeElmLtVtxEdgeSelector(mesh, wtag, maxVtx, maxEdge);

        double avgSides = parma::avgSharedSides(s);
        monitorUpdate(maxElmImb, iS, iA);
        monitorUpdate(avgSides, sS, sA);
        if( !PCU_Comm_Self() && verbose )
          fprintf(stdout, "elmImb %f avgSides %f\n", maxElmImb, avgSides);
        parma::BalOrStall* stopper =
          new parma::BalOrStall(iA, sA, sideTol*.001, verbose);

        parma::Stepper b(mesh, factor, s, w[2], t, sel, stopper);
        return b.step(tolerance, verbose);
      }
  };
}

namespace parma {
  apf::Balancer* makeElmLtVtxEdgeBalancer(apf::Mesh* m, double maxVtx,
      double maxEdge, double stepFactor, int verbosity) {
    return new ElmLtVtxEdge(m, stepFactor, maxVtx, maxEdge, verbosity);
  }
}

#include <PCU.h>
#include <stdio.h>
#include <limits.h>
#include "parma.h"
#include "parma_balancer.h"
#include "parma_step.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"
#include "parma_commons.h"
#include "parma_convert.h"

namespace {
  using parmaCommons::status;

  int getMaxNb(parma::Sides* s) {
    return PCU_Max_Int(s->size());
  }

  class ImbOrMaxNeighbor : public parma::Stop {
    public:
      ImbOrMaxNeighbor(parma::Average* nbAvg, double maxNbTol, int targets, int v=0)
        : nb(nbAvg), nbTol(maxNbTol), tgts(targets), verbose(v) {}
      bool stop(double imb, double maxImb) {
        int maxTgts = PCU_Max_Int(tgts);
        const double nbSlope = nb->avg();
        if( !PCU_Comm_Self() && verbose )
          status("max neighbor slope %f tolerance %f\n", nbSlope, nbTol);
        if( !PCU_Comm_Self() && verbose && maxTgts == 0 )
          status("no targets found... stopping\n");
        return imb > maxImb || ( fabs(nbSlope) < nbTol ) || ( maxTgts == 0 );
      }
    private:
      parma::Average* nb;
      double nbTol;
      int tgts;
      int verbose;
  };
  
  class ShapeOptimizer : public parma::Balancer {
    public:
      ShapeOptimizer(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "gap") {
      }

      bool runStep(apf::MeshTag* wtag, double tolerance) {
        parma::Sides* s = parma::makeVtxSides(mesh);
        parma::Weights* w =
          parma::makeEntWeights(mesh, wtag, s, mesh->getDimension());
        parma::Targets* t =
          parma::makeShapeTargets(s);
        parma::Selector* sel = parma::makeShapeSelector(mesh, wtag);
        double maxNb = TO_DOUBLE(getMaxNb(s));
        monitorUpdate(maxNb, sS, sA);
        parma::Stop* stopper =
          new ImbOrMaxNeighbor(sA, maxNb*.001, t->size(), verbose);
        parma::Stepper b(mesh, factor, s, w, t, sel, "elm", stopper);
        return b.step(tolerance, verbose);
      }
  };
}

apf::Balancer* Parma_MakeShapeOptimizer(apf::Mesh* m,
    double stepFactor, int verbosity) {
  return new ShapeOptimizer(m, stepFactor, verbosity);
}

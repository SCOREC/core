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

  class ImbOrLong : public parma::Stop {
    public:
      ImbOrLong(apf::Mesh* m, int tol)
        : mesh(m), sideTol(tol) {}
      bool stop(double imb, double maxImb) {
        const int small = Parma_GetSmallestSideMaxNeighborParts(mesh);
        if (!PCU_Comm_Self())
          status("Smallest Side %d, Target Side %f\n", small, sideTol);
        return imb > maxImb || small >= sideTol;
      }
    private:
      apf::Mesh* mesh;
      double sideTol;
  };
  
  class ShapeOptimizer : public parma::Balancer {
    public:
      ShapeOptimizer(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "gap") {
        smallestTgtSide = 10;
        iter=0;
        if (!PCU_Comm_Self())
          status("Factor %f Smallest target side %d\n",f, smallestTgtSide);
      }

      bool runStep(apf::MeshTag* wtag, double tolerance) {
        if (!PCU_Comm_Self())
          status("Iteration: %d\n",iter);
        parma::Sides* s = parma::makeVtxSides(mesh);
        parma::Weights* w =
          parma::makeEntWeights(mesh, wtag, s, mesh->getDimension());
        parma::Targets* t =
          parma::makeShapeTargets(s);
        parma::Selector* sel = parma::makeShapeSelector(mesh, wtag);
        ImbOrLong* stopper = new ImbOrLong(mesh, smallestTgtSide);
        parma::Stepper b(mesh, factor, s, w, t, sel, "elm", stopper);
        return b.step(tolerance, verbose);
      }
    private:
      int smallestTgtSide;
      int iter;
      int misNumber;
      int maxMis;
  };
}

apf::Balancer* Parma_MakeShapeOptimizer(apf::Mesh* m,
    double stepFactor, int verbosity) {
  return new ShapeOptimizer(m, stepFactor, verbosity);
}

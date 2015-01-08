#include <PCU.h>
#include "parma.h"
#include "parma_step.h"
#include "parma_balancer.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"

namespace {
  class VtxBalancer : public parma::Balancer {
    private:
      int sideTol;
      int getSideTol(parma::Sides* s) {
        double avg = static_cast<double>(s->total());
        PCU_Add_Doubles(&avg, 1);
        avg /= PCU_Comm_Peers();
        int tol = static_cast<int>(avg);
        if( !PCU_Comm_Self() )
          fprintf(stdout, "sideTol %d\n", tol);
        return tol;
      }
    public:
      VtxBalancer(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "vertices") {
          parma::Sides* s = parma::makeVtxSides(mesh);
          sideTol = getSideTol(s);
      }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        parma::Sides* s = parma::makeVtxSides(mesh);
        parma::Weights* w = parma::makeEntWeights(mesh, wtag, s, 0);
        parma::Targets* t =
          parma::makeWeightSideTargets(s, w, sideTol, factor);
        parma::Selector* sel = parma::makeVtxSelector(mesh, wtag);
        parma::Stepper b(mesh, factor, s, w, t, sel);
        return b.step(tolerance, verbose);
      }
  };
}

apf::Balancer* Parma_MakeVtxBalancer(apf::Mesh* m,
    double stepFactor, int verbosity) {
  if( !PCU_Comm_Self() && verbosity )
    fprintf(stdout,"PARMA_STATUS stepFactor %.3f\n", stepFactor);
  return new VtxBalancer(m, stepFactor, verbosity);
}

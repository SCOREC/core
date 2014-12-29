#include <apfPartition.h>
#include <parma.h>
#include "parma_balancer.h"
#include "parma_sides.h"
#include "parma_surfToVol.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"
#include "parma_step.h"

namespace {
  class ElmLtVtx : public parma::Balancer {
    private: 
      double maxVtx;
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
      ElmLtVtx(apf::Mesh* m, double f, double maxV, int v)
        : Balancer(m, f, v, "elements") {
          maxVtx = maxV;
          if( !PCU_Comm_Self() ) {
            fprintf(stdout, "PARMA_STATUS stepFactor %.3f\n", f);
            fprintf(stdout, "PARMA_STATUS maxVtx %.3f\n", maxVtx);
          }
        }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        const double maxVtxImb =
          Parma_GetWeightedEntImbalance(mesh, wtag, 0);
        if( !PCU_Comm_Self() )
          fprintf(stdout, "vtx imbalance %.3f\n", maxVtxImb);
        parma::Sides* s = parma::makeVtxSides(mesh);
        static int sideTol = getSideTol(s);
        PCU_Debug_Print("%s\n", s->print("sides").c_str());
        parma::Weights* w[2] =
        {parma::makeEntWeights(mesh, wtag, s, 0),
          parma::makeEntWeights(mesh, wtag, s, mesh->getDimension())};
        parma::Targets* t =
          parma::makeVtxElmTargets(s, w, sideTol, maxVtx, factor);
        delete w[0];
        PCU_Debug_Print("%s\n", t->print("targets").c_str());
        parma::Selector* sel = parma::makeElmLtVtxSelector(mesh, wtag, maxVtx);
        parma::Stepper b(mesh, wtag, factor, s, w[1], t, sel);
        return b.step(tolerance, verbose);
      }
  };
}

class VtxElmBalancer : public parma::Balancer {
  public:
    VtxElmBalancer(apf::Mesh* m, double f, int v)
      : Balancer(m, f, v, "cake") { }
    bool runStep(apf::MeshTag*, double) { return true; }
    void balance(apf::MeshTag* wtag, double tolerance) {
      apf::Balancer* b = Parma_MakeVtxBalancer(mesh, factor, verbose);
      b->balance(wtag, tolerance);
      delete b;
      Parma_PrintPtnStats(mesh, "post vertices");
      double maxVtxW = parma::getMaxWeight(mesh, wtag, 0);
      b = new ElmLtVtx(mesh, factor, maxVtxW, verbose);
      b->balance(wtag, tolerance);
      delete b;
    }
};

apf::Balancer* Parma_MakeVtxElmBalancer(apf::Mesh* m,
    double stepFactor, int verbosity) {
  return new VtxElmBalancer(m, stepFactor, verbosity);
}

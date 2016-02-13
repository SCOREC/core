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
        double avgSides = parma::avgSharedSides(s);
        if( !PCU_Comm_Self() && verbose )
          fprintf(stdout, "avgSides %f\n", avgSides);

        parma::GhostWeights* gw =
          parma::makeGhostWeights(mesh, wtag, s, layers);
        parma::Weights* elmW = convertGhostToEntWeight(gw,3);
        parma::Weights* vtxW = convertGhostToEntWeight(gw,0);
        destroyGhostWeights(gw);

        const double vtxImb = parma::getImbalance(vtxW);
        if( !PCU_Comm_Self() && verbose )
          fprintf(stdout, "vtx imbalance %.3f\n", vtxImb);

        const double elmImb = parma::getImbalance(elmW);
        monitorUpdate(elmImb, iS, iA);
        monitorUpdate(avgSides, sS, sA);

        parma::Targets* t = parma::makeTargets(s, elmW, factor);
        parma::Selector* sel = parma::makeElmSelector(mesh, wtag);
        parma::BalOrStall* stopper =
          new parma::BalOrStall(iA, sA, sideTol*.001, verbose);
        parma::Stepper b(mesh, factor, s, elmW, t, sel, stopper);
        bool ret = b.step(tolerance, verbose);
        return ret;
      }
      int layers;
  };

  class GhostVtxLtElmBalancer : public parma::Balancer {
    private:
      int sideTol;
      int layers;
      double maxElmImb;
      int stepNum;
    public:
      GhostVtxLtElmBalancer(apf::Mesh* m, double f, int v, int l)
        : Balancer(m, f, v, "ghostVtxLtElms"), layers(l), stepNum(0)
      {
        parma::Sides* s = parma::makeVtxSides(mesh);
        sideTol = static_cast<int>(parma::avgSharedSides(s));
        delete s;
        if( !PCU_Comm_Self() && verbose )
          fprintf(stdout, "sideTol %d\n", sideTol);
      }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        parma::Sides* s = parma::makeVtxSides(mesh);
        double avgSides = parma::avgSharedSides(s);
        if( !PCU_Comm_Self() && verbose )
          fprintf(stdout, "avgSides %f\n", avgSides);

        parma::GhostWeights* gw =
          parma::makeGhostWeights(mesh, wtag, s, layers);
        parma::Weights* vtxW = convertGhostToEntWeight(gw,0);
        parma::Weights* edgeW = convertGhostToEntWeight(gw,1);
        parma::Weights* elmW = convertGhostToEntWeight(gw,mesh->getDimension());
        destroyGhostWeights(gw);

        const double elmImb = parma::getImbalance(elmW);
        const double edgeImb = parma::getImbalance(edgeW);
        if( !PCU_Comm_Self() && verbose ) {
          fprintf(stdout, "elm imbalance %.3f\n", elmImb);
          fprintf(stdout, "edge imbalance %.3f\n", edgeImb);
        }
        if( !stepNum ) //FIXME need to set the imbalance at the beginning for the primary entity
          maxElmImb = elmImb;

        monitorUpdate(elmImb, iS, iA);
        monitorUpdate(avgSides, sS, sA);
        parma::Weights* w[2] = {elmW,vtxW};
        parma::Targets* t =
          parma::makeVtxElmTargets(s, w, sideTol, maxElmImb, factor); //FIXME
        delete w[0];
        parma::Selector* sel =
          parma::makeVtxLtElmSelector(mesh, wtag, maxElmImb);
        parma::BalOrStall* stopper =
          new parma::BalOrStall(iA, sA, sideTol*.001, verbose);
        parma::Stepper b(mesh, factor, s, vtxW, t, sel, stopper);
        bool ret = b.step(tolerance, verbose);
        stepNum++;
        return ret;
      }
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
      Parma_PrintPtnStats(mesh, "post-elements", (verbose>2));
      b = new GhostVtxLtElmBalancer(mesh, factor, verbose, layers);
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

#include <apf_partition.h>
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"

namespace parma {

  class MergeTargets : public Targets {
    public:
      MergeTargets(Sides* s, Weights* w, double maxImb) 
        : Targets(s,w,0.1) 
      {
        //compute average part weight and maximum imbalance 
        //if (weight < avgWeight * maxImb && weight > 0) then
        //  run knapsack and fill in the net 
        //  (see targets.h and associative.h for container API to use for net)
      }
    private:
      MergeTargets();
  };

  apf::Migration* selectMerges(apf::Mesh* m, MergeTargets& tgts) {
    //run MIS and getMergeTargets(...) to determine which 'target'
    // part this part will be merged into then create a Migration 
    // object and set the 'target' part id for each element
    //return the migration object
  };

  class MergeWeights : public EntWeights {
    public:
      MergeWeights(apf::Mesh* m, apf::MeshTag* w, Sides* s, apf::Migration* p)
        : EntWeights(m, w, s, m->getDimension()), plan(p) {}
    private:
      apf::Migration* plan;
      MergeWeights();
      double getEntWeight(apf::Mesh* m, apf::MeshEntity* e, apf::MeshTag* w) {
        assert(m->hasTag(e,w));
        double w = 0;
        if (!plan->has(e))
          m->getDoubleTag(e,w,&w);
        return w;
      }
  };

  bool canSplit(Weights& weights, double tgtImb) {
    //compute the number of empty parts
    //compute the number of empty parts needed to split parts with an 
    // imbalance greater than tgtImb s.t. their imbalance is less than
    // tgtImb
    //if ( number of empty parts >= number of parts neeeded ) 
    //  return true
    //else 
    //  return false
  }

  double imbalance(Weights* w) {
    double maxImb = 0;
    // maxImb = max weight / averge weight
    return maxImb;
  }

  double chi(apf::Mesh* m, Sides* s, Weights* w) {
    double testImb = imbalance(w);
    double step = 0.2;
    bool splits = false;
    do {
      testImb -= step;
      MergeTargets mergeTargets(sides, weights, testImb);
      Migration* plan = selectMerges(m, mergeTargets);
      MergeWeights mergeWeights(m, w, sides, plan);
      splits = canSplit(mergeWeights, testImb);
      delete plan; // not migrating
    } while ( splits );
    return testImb;
  }

  void split(Weights& weights, double tgtImb, apf::Migration* plan) {
    const int partId = PCU_Comm_Self();
    int type = 0; // { -1 if empty part, 1 if heavy part, 0 o.w.}
    int needed = 0;// p[1] = number of empty parts needed by this (heavy) part to reach tgtImb
    int payload[3] = {type, needed, partId};
    // MPI_Exscan with custom type and operator
    // see http://sbel.wisc.edu/Courses/ME964/2012/Lectures/lecture0419.pdf
    // custom operator http://www.mcs.anl.gov/research/projects/mpi/mpi-standard/mpi-report-1.1/node80.htm
    // ---write pseudo code for operator----
    // ---add elements to plan if assignment---
  }

  void hps(apf::Mesh* m, Sides* s, Weights* w, double imb) {
    MergeTargets mergeTargets(sides, weights, imb);
    Migration* plan = selectMerges(m, mergeTargets);
    MergeWeights mergeWeights(m, w, sides, plan);
    split(mergeWeights, imb, plan);
    m->migrate(plan);
  }

  class HpsBalancer : public apf::Balancer {
    public:
      HpsBalancer(apf::Mesh* m, double f, int v)
        : mesh(m), factor(f), verbose(v) 
      {
        (void) verbose; // silence!
      }
      void run(apf::MeshTag* weights, double tolerance) {
        Sides* sides = makeElmBdrySides(m);
        Weights* weights = makeEntWeights(m, w, sides, m->getDimension());
        double imb = chi(m, sides, weights);
        hps(m, sides, weights, imb);
      }
      virtual void balance(apf::MeshTag* weights, double tolerance) {
        double t0 = MPI_Wtime();
        run(weights,tolerance);
        double t1 = MPI_Wtime();
        if (!PCU_Comm_Self())
          printf("elements balanced to %f in %f seconds\n", tolerance, t1-t0);
      }
    private:
      apf::Mesh* mesh;
      double factor;
      int verbose;
  };
}; //end parma namespace

apf::Balancer* Parma_MakeHpsBalancer(apf::Mesh* m, 
    double stepFactor, int verbosity) {
  return new parma::HpsBalancer(m, stepFactor, verbosity);
}

#include <apfPartition.h>
#include <PCU.h>
#include "parma_sides.h"
#include "parma_entWeights.h"
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
      double total() {
        //return the total number of targets
        return 0;
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
        double entW = 0;
        if (!plan->has(e))
          m->getDoubleTag(e,w,&entW);
        return entW;
      }
  };

  int numSplits(Weights& w, double tgtWeight) {
    return static_cast<int>(ceil(w.self()/tgtWeight));
  }

  int isEmpty(Weights& w) {
    return (w.self() == 0) ? 1 : 0; //FIXME - dangerous comparison
  }

  int numHeavy(Weights& w, double tgtWeight) {
    int splits = numSplits(w, tgtWeight);
    PCU_Add_Ints(&splits, 1);
    return splits;
  }

  int numEmpty(Weights& w) {
    int empty = isEmpty(w);
    PCU_Add_Ints(&empty, 1);
    return empty;
  }

  bool canSplit(Weights& w, double tgt, int& extra) {
    extra = numHeavy(w, tgt) - numEmpty(w);
    if ( extra >= 0 )
      return true;
    else 
      return false;
  }

  double avgWeight(Weights* w) {
    double avg = w->self();
    PCU_Add_Doubles(&avg, 1);
    return avg/PCU_Comm_Peers();
  }

  double maxWeight(Weights* w) {
    double max = w->self();
    PCU_Max_Doubles(&max, 1);
    return max;
  }

  double imbalance(Weights* w) {
    return maxWeight(w)/avgWeight(w);
  }

  double chi(apf::Mesh* m, apf::MeshTag* wtag, Sides* s, Weights* w) {
    double testW = imbalance(w)*avgWeight(w);
    double step = 0.2;
    bool splits = false;
    int extraEmpties = 0;
    do {
      testW -= step;
      MergeTargets mergeTgts(s, w, testW);
      apf::Migration* plan = selectMerges(m, mergeTgts);
      MergeWeights mergeWeights(m, wtag, s, plan);
      splits = canSplit(mergeWeights, testW, extraEmpties);
      delete plan; // not migrating
    } while ( splits );
    assert(0==extraEmpties);
    return testW;
  }

  void split(Weights& w, double tgt, apf::Migration* plan) {
    const int partId = PCU_Comm_Self();
    int numSplit = numSplits(w, tgt);
    int empty = isEmpty(w);
    assert(!(numSplit && empty));
    int hl[2] = {numSplit, empty};
    //number the heavies and empties
    PCU_Exscan_Ints(hl, 2);
    //send heavy part ids to brokers
    PCU_Comm_Begin();
    for(int i=0; i<numSplit; i++)
      PCU_COMM_PACK(hl[0]+i, partId);
    PCU_Comm_Send();
    int heavyPartId = 0;
    int count = 0;
    while(PCU_Comm_Listen()) {
      count++;
      PCU_COMM_UNPACK(heavyPartId);
    }
    assert(count==1);
    //send empty part ids to brokers
    PCU_Comm_Begin();
    if ( empty )
      PCU_COMM_PACK(hl[1], partId); 
    PCU_Comm_Send();
    int emptyPartId = -1;
    count = 0;
    while(PCU_Comm_Listen()) {
      count++;
      PCU_COMM_UNPACK(emptyPartId);
    }
    assert(count==1);
    //brokers send empty part assignment to heavies
    PCU_Comm_Begin();
    if ( emptyPartId != -1 )
      PCU_COMM_PACK(heavyPartId, emptyPartId); 
    PCU_Comm_Send();
    std::vector<int> tgtEmpties;
    while(PCU_Comm_Listen()) {
      int tgtPartId = 0;
      PCU_COMM_UNPACK(emptyPartId);
      tgtEmpties.push_back(tgtPartId);
    }
    assert( numSplit && tgtEmpties.size() );
    //run async rib 
    //assign rib blocks to tgtEmpties 
  }

  void hps(apf::Mesh* m, apf::MeshTag* wtag, Sides* s, Weights* w, double tgt) {
    MergeTargets mergeTargets(s, w, tgt);
    apf::Migration* plan = selectMerges(m, mergeTargets);
    MergeWeights mergeWeights(m, wtag, s, plan);
    split(mergeWeights, tgt, plan);
    m->migrate(plan);
  }

  class HpsBalancer : public apf::Balancer {
    public:
      HpsBalancer(apf::Mesh* m, double f, int v)
        : mesh(m), factor(f), verbose(v) 
      {
        (void) verbose; // silence!
      }
      void run(apf::MeshTag* wtag, double tolerance) {
        Sides* sides = makeElmBdrySides(mesh);
        Weights* w = makeEntWeights(mesh, wtag, sides, mesh->getDimension());
        double tgt = chi(mesh, wtag, sides, w);
        hps(mesh, wtag, sides, w, tgt);
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

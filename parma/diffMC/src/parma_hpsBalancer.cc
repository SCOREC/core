#include <apfPartition.h>
#include <PCU.h>
#include "parma_sides.h"
#include "parma_entWeights.h"
#include "parma_targets.h"
#include "parma_selector.h"

namespace parma {

  class MergeTargets : public Targets {  // we don't really need a map/associative container here - a list/vector/array would work
    public:
      MergeTargets(Sides* s, Weights* w, double maxW) 
        : Targets(s,w,0.1) 
      {
        //if (weight < maxW && weight > 0) then
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
    return new apf::Migration(m);
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
      //getMergedWeight() // how much weight is being merged into myself??
  };

  int splits(Weights& w, double tgtWeight) {
    return static_cast<int>(ceil(w.self()/tgtWeight))-1; //FIXME - self needs to return merged part size
  }

  int isEmpty(Weights& w) {
    return (w.self() == 0) ? 1 : 0; //FIXME - dangerous comparison
  }

  int totSplits(Weights& w, double tgtWeight) {
    int numSplits = splits(w, tgtWeight);
    PCU_Add_Ints(&numSplits, 1);  // MPI_All_reduce(...,MPI_SUM,...)
    return numSplits;
  }

  int numEmpty(Weights& w) {
    int empty = isEmpty(w);
    PCU_Add_Ints(&empty, 1);
    return empty;
  }

  bool canSplit(Weights& w, double tgt, int& extra) {
    extra = numEmpty(w) - totSplits(w, tgt);
    if ( extra < 0 )
      return false;   
    else
      return true;
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

  double imbalance(apf::Mesh* m, apf::MeshTag* wtag) {
    Sides* s = makeElmBdrySides(m);
    Weights* w = makeEntWeights(m, wtag, s, m->getDimension());
    double imb = imbalance(w);
    delete w; 
    delete s;
    return imb;
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
      MergeWeights mergeWeights(m, wtag, s, plan); // compute the weight of each part post merges
      splits = canSplit(mergeWeights, testW, extraEmpties);
      delete plan; // not migrating
    } while ( splits );
    assert(0==extraEmpties);
    return testW;
  }

  void split(Weights& w, double tgt, apf::Migration* plan) {
    const int partId = PCU_Comm_Self();
    int numSplit = splits(w, tgt);
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
    //TODO run async rib 
    //TODO assign rib blocks/sub-parts to tgtEmpties
    //TODO add element empty assignments to plan 
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
      HpsBalancer(apf::Mesh* m, int v)
        : mesh(m), verbose(v) 
      {
        (void) verbose; // silence!
      }
      void run(apf::MeshTag* wtag) {
        Sides* sides = makeElmBdrySides(mesh);
        Weights* w = makeEntWeights(mesh, wtag, sides, mesh->getDimension());
        double tgt = chi(mesh, wtag, sides, w);
        hps(mesh, wtag, sides, w, tgt);
        delete sides;
        delete w;
      }
      virtual void balance(apf::MeshTag* weights, double tolerance) {
        (void) tolerance; // shhh
        double t0 = MPI_Wtime();
        run(weights);
        double elapsed = MPI_Wtime()-t0;
        PCU_Max_Doubles(&elapsed, 1);
        double maxImb = imbalance(mesh, weights);
        if (!PCU_Comm_Self())
          printf("elements balanced to %f in %f seconds\n", maxImb, elapsed);
      }
    private:
      apf::Mesh* mesh;
      int verbose;
  };
}; //end parma namespace

apf::Balancer* Parma_MakeHpsBalancer(apf::Mesh* m, int verbosity) {
  return new parma::HpsBalancer(m, verbosity);
}

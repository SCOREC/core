#include <PCU.h>
#include <apfPartition.h>
#include "parma.h"
#include "parma_sides.h"
#include "parma_entWeights.h"
#include "parma_targets.h"
#include "parma_selector.h"
#include "zeroOneKnapsack.h"
#include "maximalIndependentSet/mis.h"
#include "../zoltan/apfZoltan.h"
#include <limits>

using std::vector;
using std::set;
using parma::Weights;
using parma::Sides;

namespace {

class MergeTargets {
  public:
    //maxW == HeavyImb
    //Produces the optimal merging combination of the part's neighbors into
    // itself for a given maxW
    MergeTargets(Sides* s, Weights* w, double maxW) {
      if (w->self() >= maxW || w->self() == 0)
        return;

      //FIXME move this and soln conversion to function {
      //Create array of neighbor partId's for assignment in mergingNet
      int* nborPartIds = new int[w->size()];
      const Weights::Item* weight;
      unsigned int weightIdx = 0;
      w->begin();
      while( (weight = w->iterate()) )
        nborPartIds[weightIdx++] = weight->first;
      w->end();
      // end FIXME }

      double minWeight;
      size_t* normWeights = normalizeWeights(w, minWeight);

      const double wcap = (maxW - w->self());
      const double kdcap = floor(scale(wcap/minWeight));
      const size_t kcap = static_cast<size_t>(kdcap);


      size_t* value = new size_t[s->total()];
      std::fill(value, value + s->total(),1);

      Knapsack k = makeKnapsack(kcap, w->size(), normWeights, value);
      size_t solnVal = solve(k);
      size_t solnSz;
      size_t* soln = getSolution(k, &solnSz);
      assert(solnVal == solnSz);
      destroyKnapsack(k);

      //FIXME name too long
      mergeTargetsResults.reserve(solnSz);
      for(size_t i=0; i<solnSz; i++)
        mergeTargetsResults.push_back(nborPartIds[soln[i]]);

      free(soln);
      delete [] nborPartIds;
      delete [] value;
      delete [] normWeights;
    }

    size_t total() {
      return mergeTargetsResults.size();
    }

    int mergeTargetIndex(size_t index) {
      return mergeTargetsResults.at(index);
    }

  private:
    vector<int> mergeTargetsResults;

    double scale(double v) {
      double divideFactor = .1;
      return v/divideFactor;
    }

    size_t* normalizeWeights(Weights* w, double& minWeight){
      minWeight = std::numeric_limits<double>::max();
      const Weights::Item* weight;
      w->begin();
      while( (weight = w->iterate()) )
        if ( weight->second < minWeight )
          minWeight = weight->second;
      w->end();

      //normalizing the neighbor weights to the 
      // minimum neighbor weight with a
      // dividing factor to increase the accuracy of knapsack
      size_t* normWeights = new size_t[w->size()];
      int weightIdx = 0;
      w->begin();
      while( (weight = w->iterate()) ){
        double normW = weight->second / minWeight;
        normW = ceil(scale(normW));
        normWeights[weightIdx++] = static_cast<size_t>(normW);
      }
      w->end();

      return normWeights;
    }
};

  void generatemMisPart(apf::Mesh* m, Sides* s, MergeTargets& tgts,
    misLuby::partInfo& part){
      //Generating misLuby part info for current part
      part.id = m->getId();

      //Passing in the adjPartIds
      const Sides::Item* partId;
      s->begin();
      while( (partId = s->iterate()) )
        part.adjPartIds.push_back(partId->first);
      s->end();

      //Passing in the mergingNet
      for(size_t i = 0; i < tgts.total(); ++i)
        part.net.push_back(tgts.mergeTargetIndex(i));
      part.net.push_back(part.id);
  }

  //Using MIS Luby, selects the merges to be executed based on the merging nets
  // and returns it in the plan
  apf::Migration* selectMerges(apf::Mesh* m, Sides* s, MergeTargets& tgts) {
    misLuby::partInfo part;

    generatemMisPart(m,s,tgts,part);

    unsigned int seed = static_cast<unsigned int>(PCU_Comm_Self()+1);
    mis_init(seed);
    const int isInMis = mis(part);


    PCU_Comm_Begin();
    //If the current part is in the MIS, send notification to its mergingNet
    if (isInMis) {
      for(size_t i=0; i < tgts.total(); ++i){
        PCU_Comm_Pack(tgts.mergeTargetIndex(i), NULL, 0);
      }
    }

    PCU_Comm_Send();
    //Listening to see if it needs to merge into another part, one msg only
    //Only data needed is sender
    bool received = false;
    int destination = -1;
    while (PCU_Comm_Listen()) {
      assert(!received);
      destination = PCU_Comm_Sender();
      received = true;
    }

    //Migration of part entities to its destination part if received one
    apf::MeshIterator* it = m->begin(m->getDimension());
    apf::MeshEntity* e;
    apf::Migration* plan = new apf::Migration(m);
    while ((e = m->iterate(it)) && received)
      plan->send(e, destination);
    m->end(it);

    return plan;
  }


  int splits(Weights* w, double tgtWeight) {
    return static_cast<int>(ceil(w->self()/tgtWeight))-1;
  }

  int isEmpty(apf::Mesh* m, apf::Migration* plan) {
    const size_t numElms = m->count(m->getDimension());
    const size_t planElms = static_cast<size_t>(plan->count());
    return ((numElms - planElms) == 0) ? 1 : 0;
  }

  int numSplits(Weights* w, double tgtWeight) {
    int ns = splits(w, tgtWeight);
    PCU_Add_Ints(&ns, 1);
    return ns;
  }

  int numEmpty(apf::Mesh* m, apf::Migration* plan) {
    int empty = isEmpty(m, plan);
    PCU_Add_Ints(&empty, 1);

    return empty;
  }

  bool canSplit(apf::Mesh* m, Weights* w, apf::Migration* plan, double tgt) {
    int verbose = 1;
    const int empties = numEmpty(m, plan);
    const int splits = numSplits(w, tgt);
    const int extra = empties - splits;
    if( !PCU_Comm_Self() && verbose )
      fprintf(stdout, "HPS_STATUS chi <imb empties splits> %.3f %6d %6d\n",
        tgt, empties, splits);
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
    Sides* s = parma::makeElmBdrySides(m);
    Weights* w = parma::makeEntWeights(m, wtag, s, m->getDimension());
    double imb = imbalance(w);
    delete w;
    delete s;
    return imb;
  }
  //Function used to find the optimal heavy imbalance for executing HPS
  double chi(apf::Mesh* m, apf::MeshTag*, Sides* s, Weights* w) {
    double t0 = PCU_Time();
    int verbose = 1;
    double testW = maxWeight(w);
    double step = 0.05 * avgWeight(w);
    bool splits = false;
    do {
      testW -= step;
      MergeTargets mergeTgts(s, w, testW);
      apf::Migration* plan = selectMerges(m, s, mergeTgts);
      splits = canSplit(m, w, plan, testW);
      delete plan; // not migrating
   } while ( splits );
   testW += step;
   if( !PCU_Comm_Self() && verbose )
     fprintf(stdout, "HPS_STATUS chi ran in seconds %.3f\n", PCU_Time()-t0);
   return testW;
  }

  apf::Migration* splitPart(apf::Mesh* m, apf::MeshTag* w, int factor) {
    apf::Migration* splitPlan = NULL;
    if( !factor ) return splitPlan;
    bool isSync = false; bool isDebug = false;
    apf::Splitter* s =
      makeZoltanSplitter(m, apf::GRAPH, apf::PARTITION, isDebug, isSync);
    double imb = 1.05;
    splitPlan = s->split(w, imb, factor+1);
    delete s;
    return splitPlan;
  }


  void assignSplits(std::vector<int>& tgts,
      apf::Migration* plan) {
    assert( plan->count() );
    for (int i = 0; i < plan->count(); ++i) {
      apf::MeshEntity* e = plan->get(i);
      size_t p = static_cast<size_t>(plan->sending(e));
      assert(p <= tgts.size());
      plan->send(e, tgts[p-1]);
    }
  }

  int splits(apf::Mesh* m, Weights* w, double tgtWeight, apf::Migration* plan) {
    int extra = numEmpty(m, plan) - numSplits(w, tgtWeight);
    int split = splits(w, tgtWeight);
    int empty = isEmpty(m, plan);
    while (extra != 0) {
      double maxW = w->self();
      if( split || empty )
        maxW = 0;
      PCU_Max_Doubles(&maxW, 1);
      int id = m->getId();
      if( split || empty || maxW != w->self() )
        id = -1;
      PCU_Max_Ints(&id, 1);
      if( m->getId() == id ) {
        split = 1;
      }
      extra--;
    }
    return split;
  }

  void writePlan(apf::Migration* plan) {
    typedef std::map<int,int> mii;
    mii dest;
    for (int i = 0; i < plan->count(); ++i) {
      apf::MeshEntity* e = plan->get(i);
      dest[plan->sending(e)]++;
    }
  }

  void split(apf::Mesh* m, apf::MeshTag* wtag, Weights* w, double tgt,
      apf::Migration** plan) {
    int verbose = 1;
    assert(*plan);
    const int partId = m->getId();
    int split = splits(m, w, tgt, *plan);
    assert(split == 0 || split == 1);
    int totSplit = split;
    PCU_Add_Ints(&totSplit, 1);
    int empty = isEmpty(m, *plan);
    int totEmpty = numEmpty(m,*plan);
    if( !PCU_Comm_Self() && verbose )
      fprintf(stdout, "HPS_STATUS numEmpty %d numSplits %d\n", totEmpty, totSplit);
    assert( totSplit == totEmpty );
    assert(!(split && empty));
    int hl[2] = {split, empty};
    //number the heavies and empties
    PCU_Exscan_Ints(hl, 2);
    //send heavy part ids to brokers
    PCU_Comm_Begin();
    if( split ) {
      PCU_COMM_PACK(hl[0], partId);
    }
    PCU_Comm_Send();
    int heavyPartId = 0;
    int count = 0;
    while(PCU_Comm_Listen()) {
      count++;
      PCU_COMM_UNPACK(heavyPartId);
    }
    //send empty part ids to brokers
    PCU_Comm_Begin();
    if( empty ) {
      PCU_COMM_PACK(hl[1], partId);
    }
    PCU_Comm_Send();
    int emptyPartId = -1;
    count = 0;
    while(PCU_Comm_Listen()) {
      count++;
      PCU_COMM_UNPACK(emptyPartId);
    }
    //brokers send empty part assignment to heavies
    PCU_Comm_Begin();
    if ( emptyPartId != -1 ) {
      PCU_COMM_PACK(heavyPartId, emptyPartId);
    }
    PCU_Comm_Send();
    std::vector<int> tgtEmpties;
    while(PCU_Comm_Listen()) {
      int tgtPartId = 0;
      PCU_COMM_UNPACK(tgtPartId);
      tgtEmpties.push_back(tgtPartId);
    }
    if( split ) {
      delete *plan;
      *plan = splitPart(m, wtag, split);
      assignSplits(tgtEmpties, *plan);
    }
    writePlan(*plan);
  }

  void hps(apf::Mesh* m, apf::MeshTag* wtag, Sides* s, Weights* w, double tgt) {
    MergeTargets mergeTargets(s, w, tgt);
    apf::Migration* plan = selectMerges(m, s, mergeTargets);
    split(m, wtag, w, tgt, &plan);
    m->migrate(plan);
  }
}

namespace parma {
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
        hps(mesh, wtag, sides, w, tgt); //**Uncomment when done with testing
        delete sides;
        delete w;
        return;
      }
      virtual void balance(apf::MeshTag* weights, double tolerance) {
        (void) tolerance; // shhh
        double initImb = imbalance(mesh, weights);
        double t0 = PCU_Time();
        run(weights);
        double elapsed = PCU_Time()-t0;
        PCU_Max_Doubles(&elapsed, 1);
        double finalImb = imbalance(mesh, weights);
        if (!PCU_Comm_Self())
          printf("elements balanced from %f to %f in %f seconds\n",
              initImb, finalImb, elapsed);
      }
    private:
      apf::Mesh* mesh;
      int verbose;
  };
} //end parma namespace

apf::Balancer* Parma_MakeHpsBalancer(apf::Mesh* m, int verbosity) {
  return new parma::HpsBalancer(m, verbosity);
}

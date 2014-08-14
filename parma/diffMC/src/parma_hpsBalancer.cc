#include <apfPartition.h>
#include <PCU.h>
#include "parma_sides.h"
#include "parma_entWeights.h"
#include "parma_targets.h"
#include "parma_selector.h"
#include "zeroOneKnapsack.h"
#include "maximalIndependentSet/mis.h"
#include <limits>

using std::vector;
using std::set;

namespace parma {

void print(int level, char* printStr,...) {
  if (level == 2) {
    // PCU_Debug_Print(printStr, ...);
  }
  else if (!PCU_Comm_Self()){}
    // printf(printStr, ...);
}

class MergeTargets {
  public:
    //maxW == HeavyImb
    //Produces the optimal merging combination of the part's neighbors into
    // itself for a given maxW
    MergeTargets(Sides* s, Weights* w, double maxW, bool chi) {
      //If part is heavy or empty, exit function
      if (w->self() >= maxW || w->self() == 0) {
        //PCU_Debug_Print("Part %d of weight %f is heavy\n", PCU_Comm_Self(),
            // w->self()); rating 0 (only in HPS, not CHI)
        return;
      }
      /*PCU_Debug_Print("Part %d of weight %f is light"
          "and not empty with imb of %f, %zu neighbors--\n",
          PCU_Comm_Self(), w->self(), maxW, w->size());
      rating 0 (only in HPS, not CHI)*/

      //Create array of neighbor partId's for assignment in mergingNet
      // PCU_Debug_Print("Neighbor part Ids \n"); //rating 1
      int* nborPartIds = new int[w->size()];
      const Weights::Item* weight;
      unsigned int weightIdx = 0;
      w->begin();
      while( (weight = w->iterate()) )
        nborPartIds[weightIdx++] = weight->first;
        // PCU_Debug_Print("\tNeighbor %d = %d\n", weightIdx-1,
        //   nborPartIds[weightIdx-1]); //rating 1
      w->end();

      double minWeight = std::numeric_limits<double>::max();
      int* normalizedIntWeights = normalizeWeights(w, minWeight);

      //How much weight can be added on to the current part before it
      // reaches the HeavyImb
      const double weightCapacity = (maxW - w->self());
      //total weight that can be added to self, normalized to min
      const int knapsackCapacity = floor(scale(weightCapacity/minWeight));

      /*PCU_Debug_Print("Weight Capcity = %f \nknapsackCapacity = %d \n",
          weightCapacity, knapsackCapacity);
      rating 1*/
      //Knapsack execution
      //Declared for knapsack class
      int* value = new int[s->total()];
      std::fill (value, value + s->total(),1);

      knapsack* ks = new knapsack(knapsackCapacity, w->size(),
          normalizedIntWeights, value);
      const int solnVal = ks->solve();
      mergeTargetsResults.reserve(solnVal);
      ks->getSolution(mergeTargetsResults);

      //PCU_Debug_Print("mergetargets size = %zu\n", mergeTargetsResults.size()); rating 0 (1 if in chi)

      //Converting mergeTargetsResults to partId's
      vector<int> partIdMergeTargets;
      partIdMergeTargets.reserve(mergeTargetsResults.size());
      for(size_t i=0; i<mergeTargetsResults.size(); i++)  {
        partIdMergeTargets.push_back(nborPartIds[mergeTargetsResults[i]]);
      }
      //Constant function to change vectors
      mergeTargetsResults.swap(partIdMergeTargets);

      //Debug to see which parts contained in the results
      for(size_t i=0; i<mergeTargetsResults.size(); i++)  {
        PCU_Debug_Print("\tmerge Target result %d\n", mergeTargetsResults[i]); //rating 1 (2 if in chi)
      }

      delete [] nborPartIds;
      delete [] value;
      delete [] normalizedIntWeights;
      delete ks;

    }

    size_t total() {
      return mergeTargetsResults.size();
    }

    int mergeTargetIndex(int index) {
      return mergeTargetsResults.at(index);
    }

  private:
    MergeTargets();
    vector<int> mergeTargetsResults;

    double scale(double v) {
      double divideFactor = .1;
      return v/divideFactor;
    }

    int* normalizeWeights(Weights* w, double& minWeight){
      //iterating through the neighbor weights to determine the minimum weight
      const Weights::Item* weight;
      w->begin();
      while( (weight = w->iterate()) )
        if ( weight->second < minWeight )
          minWeight = weight->second;
      w->end();

      // PCU_Debug_Print("min weight == %f\n", minWeight); rating 2

      //normalizing the neighbor weights to the minimum neighbor weight with a
      // dividing factor to increase the accuracy of knapsack
      int* normalizedIntWeights = new int[w->size()];
      int weightIdx = 0;
      w->begin();
      while( (weight = w->iterate()) ){
        //Divide Factor normalizing weight code
        double normalizedWeight = weight->second / minWeight;
        normalizedWeight = scale(normalizedWeight);
        normalizedIntWeights[weightIdx++] = (int) ceil(normalizedWeight);

        // PCU_Debug_Print("weight %d, normalized to %d from %f\n",
        // weight->first, normalizedIntWeights[(weightIdx-1)], weight->second); //rating 1, 2 if CHI
      }
      w->end();

      return normalizedIntWeights;
    }
};

  void generatemMisPart(apf::Mesh* m, Sides* s, MergeTargets& tgts,
    vector<misLuby::partInfo>& parts){
      //Generating misLuby part info for current part
      misLuby::partInfo part;
      part.id = m->getId();

      //Passing in the adjPartIds
      const Sides::Item* partId;
      s->begin();
      while( (partId = s->iterate()) )
        part.adjPartIds.push_back(partId->first);
      s->end();

      PCU_Debug_Print("adjpartIds size = %zu\n", part.adjPartIds.size()); //rating 2

      PCU_Debug_Print("Part %d mergeNet size %zu\n", PCU_Comm_Self(),
        tgts.total()); //rating 0

      //Passing in the mergingNet
      for(size_t i = 0; i < tgts.total(); ++i){
        part.net.push_back(tgts.mergeTargetIndex(i));
        PCU_Debug_Print("\t%zu mergingNet %d\n",i , part.net[i]);//rating 1 (CHI 2)
      }
      part.net.push_back(part.id);

    //Testing for examples, additionally add paramter [4]
    //in mis as true if testing random numbers
    //Change random numbers as you please, MIS selects lowest numbers (only ints)
    //Additionally, more can be added for more parts
    if(part.id == 0) part.randNum = 1;
    else if (part.id == 1) part.randNum = 2;
    else if (part.id == 2) part.randNum = 3;
    else if (part.id == 3) part.randNum = 1;
    else if (part.id == 4) part.randNum = 2;
    else if (part.id == 5) part.randNum = 2;
    else if (part.id == 6) part.randNum = 1;
    else if (part.id == 7) part.randNum = 2;

      parts.push_back(part);
      //End creating misLuby part
  }

  //Using MIS Luby, selects the merges to be executed based on the merging nets
  // and returns it in the plan
  apf::Migration* selectMerges(apf::Mesh* m, Sides* s, MergeTargets& tgts) {
    //TODO Vector to be passed into misLuby (later to be changed to just
    //a single part info since vector is from multi parts per process)
    vector<misLuby::partInfo> parts;
    parts.reserve(1);

    generatemMisPart(m,s,tgts,parts);

    int randNumSeed = time(NULL)%(PCU_Comm_Self()+1);

    mis_init(randNumSeed,false);
    vector<int> maximalIndSet;
    //Maximal Independent Set Call
    int ierr = mis(PCU_Comm_Self(), s->total()+1, parts, maximalIndSet, true);
    //Assert will fail if MIS fails
    assert(!ierr);

    // Debug if in MIS rating 0
    if (maximalIndSet.size() == 1)
      PCU_Debug_Print("Part %d in MIS\n", PCU_Comm_Self()); //rank 0 (1 if CHI)
    else PCU_Debug_Print("Part %d NOT in MIS\n", PCU_Comm_Self()); //rank 0 (1 if CHI)

    PCU_Comm_Begin();
    //If the current part is in the MIS, send notification to its mergingNet
    if (maximalIndSet.size() != 0) {
      for(size_t i=0; i < tgts.total(); ++i){
        //PCU_Debug_Print("send dest = %d\n", tgts.mergeTargetIndex(i)); rating 1
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
    PCU_Debug_Print("destination = %d\n", destination); //rating 0 (1 if CHI)

    //Migration of part entities to its destination part if received one
    apf::MeshIterator* it = m->begin(m->getDimension());
    apf::MeshEntity* e;
    apf::Migration* plan = new apf::Migration(m);
    while ((e = m->iterate(it)) && received)
      plan->send(e, destination);
    m->end(it);

    return plan;
  };

  int splits(Weights* w, double tgtWeight) {
    return static_cast<int>(ceil(w->self()/tgtWeight))-1;
  }

  int isEmpty(apf::Mesh* m, apf::Migration* plan) {
    return (m->count(m->getDimension()) - plan->count() == 0) ? 1 : 0;
  }

  int totSplits(Weights* w, double tgtWeight) {
    int numSplits = splits(w, tgtWeight);
    PCU_Add_Ints(&numSplits, 1);
    return numSplits;
  }

  int numEmpty(apf::Mesh* m, apf::Migration* plan) {
    int empty = isEmpty(m, plan);
    PCU_Add_Ints(&empty, 1);

    return empty;
  }

  bool canSplit(apf::Mesh* m, Weights* w, apf::Migration* plan, double tgt,
      int& extra)
  {
    extra = numEmpty(m, plan) - totSplits(w, tgt);
    //PCU_Debug_Print("Empty parts = %d - total Splits = %d = Extra %d ",
    // numEmpty(m,plan),totSplits(w,tgt), extra); rating 2
    if ( extra < 0 ){
      return false;  }
    else {
      return true;
     }
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
  //Function used to find the optimal heavy imbalance for executing HPS
  double chi(apf::Mesh* m, apf::MeshTag* wtag, Sides* s, Weights* w) {
    double testW = maxWeight(w);
    //Step size can change arbitrarily as needed
    double step = 0.1 * avgWeight(w);
    bool splits = false;
    int extraEmpties = 0;
    do {
      testW -= step;
      //**IDEA** Add sizes of mergeNets across processes and subtract
      // from totSplits and if postitve, fail
      MergeTargets mergeTgts(s, w, testW, false);
      apf::Migration* plan = selectMerges(m, s, mergeTgts);
      splits = canSplit(m, w, plan, testW, extraEmpties);
      PCU_Debug_Print("Test imb %f, result %d\n", testW, splits); //rating 0
      delete plan; // not migrating
   }
   while ( splits );
    testW += step;
   return testW;
  }

  void split(apf::Mesh* m, Weights* w, double tgt, apf::Migration* plan) {
    const int partId = PCU_Comm_Self();
    int numSplit = splits(w, tgt);
    int empty = isEmpty(m, plan);
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
    MergeTargets mergeTargets(s, w, tgt, true);
    apf::Migration* plan = selectMerges(m, s, mergeTargets);
    split(m, w, tgt, plan);
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
        // PCU_Debug_Print("Final Chi = %f",tgt); rating 0
        //hps(mesh, wtag, sides, w, tgt); **Uncomment when done with testing
        delete sides;
        delete w;
        return;
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

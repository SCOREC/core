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
  //header for avgWeights
  double avgWeight(Weights* w);

  class MergeTargets { 
    public:
      //Have this storing the results in Targets associative class and assuming maxW = avgWeight * maxImb
      //maxW is also = HeavyImb
      MergeTargets(Sides* s, Weights* w, double maxW)
      {

        //If part is heavy or empty, exit function
        if (w->self() >= maxW || w->self() == 0) {
          PCU_Debug_Print("Part %d of weight %f is heavy\n", PCU_Comm_Self(), w->self());  
          return;
        }
        PCU_Debug_Print("Part %d of weight %f is light and not empty with imb of %f, %u neighbors--\n", PCU_Comm_Self(), w->self(), maxW, w->size());
        
        int* nborPartIds = new int[w->size()];
        int i = 0;
        const Weights::Item* weight;        
        w->begin();         
        while( (weight = w->iterate()) ) 
          nborPartIds[i++] = weight->first;    
        w->end();         

        //iterating through the neighbor weights to determine the minimum weight
        double minWeight = std::numeric_limits<double>::max(); 
        w->begin(); 
        while( (weight = w->iterate()) ) 
          if ( weight->second < minWeight )
            minWeight = weight->second;
        w->end();
        //normalizing the neighbor weights with a dividing factor to increase the accuracy of knapsack
        double t0 = MPI_Wtime();
        bool divide = true;
        int* normalizedIntWeights = new int[w->size()];
        unsigned int weightIdx = 0;
        double divide_factor = .1;
        w->begin(); 
        while( (weight = w->iterate()) ){
          
          //Divide Factor normalizing weight code
          double normalized_weight = weight->second / minWeight;
          normalized_weight /= divide_factor;
          normalizedIntWeights[weightIdx++] = (int) ceil(normalized_weight); 
          
          //interger normalizing weight code
          // normalizedIntWeights[weightIdx++] = (int) ceil( weight->second / minWeight );
          // divide = false; 

          PCU_Debug_Print("weight %d, normalized to %d from %f\n", weight->first, normalizedIntWeights[(weightIdx-1)], weight->second);
        }

        w->end();            

        const double weightCapacity = (maxW - w->self()); //How much weight can be added on to the current part before it reaches the HeavyImb
        const int knapsackCapacity = floor((weightCapacity/minWeight)/divide_factor); //total weight that can be added to self, normalized to min

        // const int knapsackCapacity = floor(weightCapacity/minWeight);

        PCU_Debug_Print("Weight Capcity = %f \nknapsackCapacity = %d \n", weightCapacity, knapsackCapacity);


        int* value = new int[s->total()];
        std::fill (value, value + s->total(),1);

        knapsack* ks = new knapsack(knapsackCapacity, w->size(), normalizedIntWeights, value);        
        const int solnVal = ks->solve();    
        mergeTargetsResults.reserve(solnVal);
        ks->getSolution(mergeTargetsResults);
        double t1 = MPI_Wtime();
        //if (divide) printf("Part %d executed knapsack in %f with knapsack type dividing factor, %f\n", PCU_Comm_Self(), t1 - t0, divide_factor);
        //else printf("Part %d executed knapsack in %f with knapsack type interger rounding\n", PCU_Comm_Self(), t1 - t0);
        
        //PCU_Debug_Print("mergetargets start\n");  
        PCU_Debug_Print("mergetargets size = %u\n", mergeTargetsResults.size());

        //Converting mergeTargetsResults to partId's
        vector<int> partIdMergeTargets;
        partIdMergeTargets.reserve(mergeTargetsResults.size());
        
        for(size_t i=0; i<mergeTargetsResults.size(); i++)  {
          partIdMergeTargets.push_back(nborPartIds[mergeTargetsResults[i]]);
          
           // PCU_Debug_Print("merge Target result %d -- ", mergeTargetsResults[i]);
           // PCU_Debug_Print("merge target part ID %d\n", nborPartIds[mergeTargetsResults[i]]);
        }

        mergeTargetsResults.swap(partIdMergeTargets);

        for(size_t i=0; i<mergeTargetsResults.size(); i++)  {
          PCU_Debug_Print("\tmerge Target result %d\n", mergeTargetsResults[i]);
        }
         
        delete [] nborPartIds;
        delete [] value;  
        delete [] normalizedIntWeights;
        delete ks;    
               
      }

      size_t total() {
        return mergeTargetsResults.size();
      }

      const int mergeTargetIndex(int& index){
        return mergeTargetsResults.at(index);
      }
      
    private:
      MergeTargets();
      vector<int> mergeTargetsResults;
  };


  //Create example of random numbers that will and will not find a merge. 2 or 3 //Done
  //1a) Look into luby paper about running maximal solution to tell you some kind of properties about the maximum //None that I could find
  //1b) Approximation of the maximum independent set using a maximal independent set undefined algo. //NPC unless graph can be classified
  //2a) Look into running it multiple times  //Far as i can tell, nothing
  //2b) **Possibly assign an upper bound to mergeTargets in Chi in that there's no possible way to get enough empties //Believe this is still doable (ask about it tomorrow)
  //3) write getMergeTargets(...) function

  //4) Create PDF with exmaples and make them highly instructional for others to run them
  apf::Migration* selectMerges(apf::Mesh* m, Sides* s, MergeTargets& tgts) {
    //run MIS and getMergeTargets(...) to determine which 'target'
    // part this part will be merged into then create a Migration 
    // object and set the 'target' part id for each element
    //return the migration object

    //Vector to be passed into misLuby (later to be changed to just a single part info since vector is still from multi parts per process)
    vector<misLuby::partInfo> parts;
    parts.reserve(1);
    
    //Generating misLuby part info for current part 
    misLuby::partInfo part;
    part.id = m->getId(); 

    //Passing in the adjPartIds (works)
    const Sides::Item* partId;  
    s->begin();         
    while( (partId = s->iterate()) ) 
      part.adjPartIds.push_back(partId->first);
    s->end();
      /*
      PCU_Debug_Print("adjPartIds:\n");
      vector<int>::iterator itr = part.adjPartIds.begin();
      while(itr != part.adjPartIds.end()){
        PCU_Debug_Print("\t%i\n", *itr++);
        }
      PCU_Debug_Print("adjPartIds end\n");
      */

    //Passing in the mergingNet (works)
    PCU_Debug_Print("Part %d mergeNet size %d\n", PCU_Comm_Self(), tgts.total());
    for(int i = 0; i < tgts.total(); ++i){
      part.net.push_back(tgts.mergeTargetIndex(i));
      PCU_Debug_Print("\tmergingNet %d\n", part.net[i]);
    }
    part.net.push_back(part.id);
    int temp_partId = part.id;
      /*
      //Testing for examples, additionally add paramter [4] in mis as true if testing random numbers
      if(temp_partId == 0) part.randNum = 1;
      else if (temp_partId == 1) part.randNum = 2;
      else if (temp_partId == 2) part.randNum = 3;
      else if (temp_partId == 3) part.randNum = 1;
      else if (temp_partId == 4) part.randNum = 2;
      else if (temp_partId == 5) part.randNum = 2;
      else if (temp_partId == 6) part.randNum = 1;
      else if (temp_partId == 7) part.randNum = 2;
      */

    parts.push_back(part); //End creating misLuby part info
    
    mis_init(0,false); //RandNumSeed set to 0 for testing purposes.
    vector<int> maximalIndSet; 
    int ierr = mis(PCU_Comm_Self(), s->total(), parts, maximalIndSet, true); //rank, total number of parts, parts mis vector, maximalindset vec
    //Todo** Remove true when done testing random numbers

    
    if (!ierr) PCU_Debug_Print("Part %d had a successful MIS\n", PCU_Comm_Self());
    if (maximalIndSet.size() == 1) PCU_Debug_Print("\tPart %d in MIS\n", part.id);
    
    PCU_Comm_Begin();

    if (maximalIndSet.size() != 0) {
      for(int i=0; i < tgts.total(); ++i){
        PCU_Debug_Print("send dest = %d\n", tgts.mergeTargetIndex(i));
        PCU_Comm_Pack(tgts.mergeTargetIndex(i), NULL, 0);
      }
    }
    
    PCU_Comm_Send();
    bool received = false;
    int destination = -1;
    while (PCU_Comm_Listen()) {
      assert(!received);
      destination = PCU_Comm_Sender();
      received = true;
    }
    PCU_Debug_Print("destination = %d\n", destination);

    apf::MeshIterator* it = m->begin(m->getDimension());
    apf::MeshEntity* e;
    apf::Migration* plan = new apf::Migration(m);
    while ((e = m->iterate(it)) && received)
      plan->send(e, destination);
    m->end(it);
    
    return plan;
    //Check to see 
  };

  int splits(Weights* w, double tgtWeight) {
    return static_cast<int>(ceil(w->self()/tgtWeight))-1; 
  }

  int isEmpty(apf::Mesh* m, apf::Migration* plan) {
    return (m->count(m->getDimension()) - plan->count() == 0) ? 1 : 0;
  }

  int totSplits(Weights* w, double tgtWeight) {
    int numSplits = splits(w, tgtWeight);
    PCU_Add_Ints(&numSplits, 1);  // MPI_All_reduce(...,MPI_SUM,...)
    return numSplits;
  }

  int numEmpty(apf::Mesh* m, apf::Migration* plan) {
    int empty = isEmpty(m, plan);
    PCU_Add_Ints(&empty, 1);

    return empty;
  }

  bool canSplit(apf::Mesh* m, Weights* w, apf::Migration* plan, double tgt, int& extra) {
    PCU_Debug_Print("Empty parts = %d, total Splits = %d $$ ", numEmpty(m,plan), totSplits(w,tgt));
    extra = numEmpty(m, plan) - totSplits(w, tgt);
    PCU_Debug_Print("Extra = %d\n", extra);
    if ( extra < 0 ){
      return false;  } 
    else{
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
  //Should we test the maxWeight first and not do a step first?
  double chi(apf::Mesh* m, apf::MeshTag* wtag, Sides* s, Weights* w) {
    double testW = maxWeight(w)*1.2;//*1.1476; 
    double step = 0.1 * avgWeight(w); //Not necessarily arbitrary but can be changed.
    bool splits = false;
    int extraEmpties = 0;
    do {
      testW -= step;
      //PCU_Debug_Print("Test Imb = %f\n", testW);
      MergeTargets mergeTgts(s, w, testW);
      // PCU_Debug_Print("Test mergeindex 0, %u",mergeTgts.mergeTargetIndex(i));
      apf::Migration* plan = selectMerges(m, s, mergeTgts); 
      splits = canSplit(m, w, plan, testW, extraEmpties);
      delete plan; // not migrating
   } while ( splits );
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
    MergeTargets mergeTargets(s, w, tgt);
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
        PCU_Debug_Print("Final Chi = %f\n======================================\n",tgt);
        delete sides;
        delete w;
        return; //TODO remove return after testing and put deletes below
        hps(mesh, wtag, sides, w, tgt);

      }
      virtual void balance(apf::MeshTag* weights, double tolerance) {
        (void) tolerance; // shhh
        double t0 = MPI_Wtime();
        PCU_Debug_Print("======================================\n");
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

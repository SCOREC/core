#include "phPartition.h"
#include "phInput.h"
#include "ph.h"
#include <parma.h>
#include <apfZoltan.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <pcu_util.h>
#include <lionPrint.h>



#ifdef HAVE_SIMMETRIX
#include <apfSIM.h>
#include "SimPartitionedMesh.h"
#include "SimModel.h"
#include "SimUtil.h"
#endif

namespace ph {

void setWeight(apf::Mesh* m, apf::MeshTag* tag, int dim) {
  double w = 1.0;
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(dim);
  int nverts =1;
  apf::Downward verts;
  while ((e = m->iterate(it))){
    int dimEnt=getDimension(m,e);
    if(dimEnt==3) {
      w=1.0;
      nverts = m->getDownward(e, 0, verts);
      if(nverts==8) w=6.0;
      if(nverts==6) w=3.0; 
    }
    m->setDoubleTag(e, tag, &w);
  }
  m->end(it);
}

apf::MeshTag* setWeights(apf::Mesh* m) {
  apf::MeshTag* tag = m->createDoubleTag("parma_weight", 1);
  setWeight(m, tag, 0);
  setWeight(m, tag, m->getDimension());
  return tag;
}


apf::Migration* getSplitPlan(Input& in, apf::Mesh2* m)
{
  PCU_ALWAYS_ASSERT(in.splitFactor >= 1);
  apf::Migration* plan;
  if (in.splitFactor != 1 ) { //0 ) { // 1) {  // 0 will force reparttition  even when split factor is 1
    apf::Splitter* splitter;
    if (in.partitionMethod == "rib") { //prefer SCOREC RIB over Zoltan RIB
      splitter = Parma_MakeRibSplitter(m);
    } else {
      std::map<std::string, int> methodMap;
      methodMap["graph"] = apf::GRAPH;
      methodMap["zrib"] = apf::RIB;
      methodMap["hypergraph"] = apf::HYPERGRAPH;
      int method = methodMap[in.partitionMethod];
      if(in.localPtn == true)
        splitter = apf::makeZoltanSplitter(m, method, apf::REPARTITION);
      else
        splitter = apf::makeZoltanGlobalSplitter(m, method, apf::REPARTITION);
    }
//    apf::MeshTag* weights = Parma_WeighByMemory(m);
    apf::MeshTag* weights = setWeights(m);
//    Parma_PrintPtnStats(m, "preSplit",false);
    plan = splitter->split(weights, 1.01, in.splitFactor);
    apf::removeTagFromDimension(m, weights, m->getDimension());
    m->destroyTag(weights);
    delete splitter;
  } else {
    plan = new apf::Migration(m);
  }
  return plan;
}

apf::Migration* split(Input& in, apf::Mesh2* m)
{
  return getSplitPlan(in,m);
}

bool isMixed(apf::Mesh2* m) {
  int mixed = 0;
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(m->getDimension());
  while ((e = m->iterate(it)))
    if ( m->getType(e) != apf::Mesh::TET ) {
      mixed = 1;
      break;
    }
  m->end(it);
  return m->getPCU()->Max<int>(mixed);
}

void clearTags(apf::Mesh* m, apf::MeshTag* t) {
  apf::removeTagFromDimension(m, t, 0);
  apf::removeTagFromDimension(m, t, m->getDimension());
}


void neighborReduction(apf::Mesh2* m, apf::MeshTag* weights, int verbose, bool fineStats) {
  double elmImb = Parma_GetWeightedEntImbalance(m, weights, m->getDimension());
  double elmImbLimit = 1+((elmImb-1)*1.20);
  double ignored = 0.1;
  apf::Balancer* balancer = Parma_MakeShapeOptimizer(m, ignored, verbose);
  balancer->balance(weights, elmImbLimit);
  delete balancer;
  Parma_PrintPtnStats(m, "postGap", fineStats);
}

void parmaMixed(Input& in, apf::Mesh2* m) {
  bool fineStats=false; // set to true for per part stats
  Parma_PrintPtnStats(m, "preRefine", fineStats); //FIXME
//  apf::MeshTag* weights = Parma_WeighByMemory(m);
  apf::MeshTag* weights = setWeights(m);
  const double step = 0.2;
  const int verbose = 0;
  apf::Balancer* balancer = Parma_MakeElmBalancer(m, step, verbose);
  balancer->balance(weights, in.elementImbalance);
  delete balancer;
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
}

void parmaTet(Input& in, apf::Mesh2* m, bool runGap) {
  bool fineStats=false; // set to true for per part stats
  Parma_PrintPtnStats(m, "preRefine", fineStats); //FIXME
  apf::MeshTag* weights = setWeights(m);
  const double step = 0.3;
  const int verbose = in.parmaVerbosity;  // set to 2 for per iteration stats
  if(runGap)
    neighborReduction(m,weights,verbose,fineStats);
  for(int i=0; i<in.parmaLoops; i++) {
    apf::Balancer* balancer = Parma_MakeVtxElmBalancer(m, step, verbose);
    balancer->balance(weights, in.vertexImbalance);
    Parma_PrintPtnStats(m, "post Parma_MakeVtxElmBalancer", fineStats);
    delete balancer;
    if(runGap)
      neighborReduction(m,weights,verbose,fineStats);
    double vtxImb = Parma_GetWeightedEntImbalance(m, weights, 0);
    if( vtxImb <= in.vertexImbalance ) {
      if( !m->getPCU()->Self() )
        lion_oprint(1, "STATUS vtx imbalance target %.3f reached\n",
            in.vertexImbalance);
      break;
    }
  }
  clearTags(m, weights);
  m->destroyTag(weights);
}

void parmaBalance(Input& in, apf::Mesh2* m, bool runGap) {
//  if( isMixed(m) )
//    parmaMixed(in,m);
//  else
    parmaTet(in,m,runGap);
}

void runBalancer(apf::Mesh2* m, Input& in, apf::Balancer* b)
{
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  b->balance(weights, in.elementImbalance);
  delete b;
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
}

void zoltanBalance(apf::Mesh2* m, Input& in, int method)
{
  bool dbg=false;
  runBalancer(m, in, apf::makeZoltanBalancer(
        m, method, apf::REPARTITION, dbg));
}

#ifdef HAVE_SIMMETRIX
void setEqualWeights(pParMesh pmesh, int desiredTotNumParts, pProgress progress)
{
  pPartitionOpts pOpts = PartitionOpts_new();
  // Set total no. of parts
  PartitionOpts_setTotalNumParts(pOpts, desiredTotNumParts);
  // Sets processes to be equally weighted
  PartitionOpts_setProcWtEqual(pOpts);
  PM_partition(pmesh, pOpts, progress);     // Do the partitioning
  PartitionOpts_delete(pOpts);              // Done with options
}

void simmetrixBalance(apf::Mesh2* m)
{
  // start progress
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  // get simmetrix mesh
  apf::MeshSIM* apf_msim = dynamic_cast<apf::MeshSIM*>(m);
  pParMesh pmesh = apf_msim->getMesh();

  // get total number of processors
  // we assume one part per processor
  int totalNumParts = PMU_size();

  // current total num parts in pmesh cannot be more than requested
  int currentTotalNumParts = PM_totalNumParts(pmesh);
  if (currentTotalNumParts > totalNumParts) {
    if( !m->getPCU()->Self() )
      lion_eprint(1, "Error: cannot reduce number of partitions %d->%d\n",
              currentTotalNumParts, totalNumParts);
    totalNumParts = currentTotalNumParts;
  }

  // set weights and do balacing
  setEqualWeights(pmesh, totalNumParts, progress);

  // delete progress
  Progress_delete(progress);
}
#endif

void balance(Input& in, apf::Mesh2* m)
{
  if(in.prePhastaBalanceMethod == "parma-gap" )
    parmaBalance(in,m,true);
  else if(in.prePhastaBalanceMethod == "parma" )
    parmaBalance(in,m,false);
  else if(in.prePhastaBalanceMethod == "zrib" )
    zoltanBalance(m,in,apf::RIB);
  else if(in.prePhastaBalanceMethod == "graph" )
    zoltanBalance(m,in,apf::GRAPH);
#ifdef HAVE_SIMMETRIX
  else if(in.prePhastaBalanceMethod == "simmetrix" )
    simmetrixBalance(m);
#endif
  else if(in.prePhastaBalanceMethod != "none" )
    fail("unknown setting for prePhastaBalanceMethod \"%s\"\n",
      in.prePhastaBalanceMethod.c_str());
}

}

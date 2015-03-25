#include <PCU.h>
#include "phPartition.h"
#include "phInput.h"
#include <parma.h>
#include <apfZoltan.h>
#include <apfMDS.h>
#include <apfMesh2.h>

namespace ph {

apf::Migration* getSplitPlan(Input& in, apf::Mesh2* m)
{
  assert(in.recursivePtn <= 1);
  assert(in.splitFactor >= 1);
  apf::Migration* plan;
  if (in.splitFactor != 1) {
    apf::Splitter* splitter;
    if (in.partitionMethod == "rib") { //prefer SCOREC RIB over Zoltan RIB
      splitter = Parma_MakeRibSplitter(m);
    } else {
      std::map<std::string, int> methodMap;
      methodMap["graph"] = apf::GRAPH;
      methodMap["hypergraph"] = apf::HYPERGRAPH;
      int method = methodMap[in.partitionMethod];
      splitter = apf::makeZoltanSplitter(m, method, apf::REPARTITION);
    }
    apf::MeshTag* weights = Parma_WeighByMemory(m);
    plan = splitter->split(weights, 1.03, in.splitFactor);
    apf::removeTagFromDimension(m, weights, m->getDimension());
    m->destroyTag(weights);
    delete splitter;
  } else {
    plan = new apf::Migration(m);
  }
  return plan;
}

void split(Input& in, apf::Mesh2* m, void (*runAfter)(apf::Mesh2*))
{
  apf::splitMdsMesh(m, getSplitPlan(in, m), in.splitFactor, runAfter);
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
  PCU_Max_Ints(&mixed, 1);
  return mixed;
}

void setWeight(apf::Mesh* m, apf::MeshTag* tag, int dim) {
  double w = 1.0;
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(dim);
  while ((e = m->iterate(it)))
    m->setDoubleTag(e, tag, &w);
  m->end(it);
}

apf::MeshTag* setWeights(apf::Mesh* m) {
  apf::MeshTag* tag = m->createDoubleTag("parma_weight", 1);
  setWeight(m, tag, 0);
  setWeight(m, tag, m->getDimension());
  return tag;
}

void clearTags(apf::Mesh* m, apf::MeshTag* t) {
  apf::removeTagFromDimension(m, t, 0);
  apf::removeTagFromDimension(m, t, m->getDimension());
}

void balance(apf::Mesh2* m)
{
  bool fineStats=false; // set to true for per part stats
  Parma_PrintPtnStats(m, "preRefine", fineStats);
  if ( isMixed(m) ) {

    apf::MeshTag* weights = Parma_WeighByMemory(m);
    double tolerance = 1.05;
    const double step = 0.2;
    const int verbose = 0;
    apf::Balancer* balancer = Parma_MakeElmBalancer(m, step, verbose);
    balancer->balance(weights, tolerance);
    delete balancer;
    apf::removeTagFromDimension(m, weights, m->getDimension());
    m->destroyTag(weights);

  } else {
    apf::MeshTag* weights = setWeights(m);
    const double vtxImbTol = 1.03;
    const double step = 0.3;
    const int verbose = 1;  // set to 2 for per iteration stats
    const double ignored = 42.42;

    Parma_ProcessDisconnectedParts(m);
    Parma_PrintPtnStats(m, "post ProcessDisconnectedParts", fineStats);

    apf::Balancer* balancer = Parma_MakeHpsBalancer(m,verbose);
    balancer->balance(weights, ignored);
    delete balancer;
    Parma_PrintPtnStats(m, "post HPS", fineStats);

    for(int i=0; i<3; i++) {
      balancer = Parma_MakeVtxElmBalancer(m, step, verbose);
      balancer->balance(weights, vtxImbTol);
      Parma_PrintPtnStats(m, "post Parma_MakeVtxElmBalancer", fineStats);
      delete balancer;
      double vtxImb = Parma_GetWeightedEntImbalance(m, weights, 0);
      if( vtxImb <= vtxImbTol ) {
        if( !PCU_Comm_Self() )
          fprintf(stdout, "STATUS vtx imbalance target %.3f reached\n",
            vtxImbTol);
        break;
      }
    }
    clearTags(m, weights);
    m->destroyTag(weights);
  }
}

}

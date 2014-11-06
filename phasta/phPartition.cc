#include <PCU.h>
#include "phPartition.h"
#include "phInput.h"
#include <parma.h>
#include <apfZoltan.h>
#include <apfMDS.h>
#include <apfMesh2.h>

namespace ph {

void split(Input& in, apf::Mesh2* m, void (*runAfter)(apf::Mesh2*))
{
  assert(in.recursivePtn <= 1);
  int factor = in.numTotParts / PCU_Comm_Peers();
  assert(in.numTotParts % PCU_Comm_Peers() == 0);
  apf::Migration* plan;
  if (factor != 1) {
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
    plan = splitter->split(weights, 1.03, factor);
    apf::removeTagFromDimension(m, weights, m->getDimension());
    m->destroyTag(weights);
    delete splitter;
  } else {
    plan = new apf::Migration(m);
  }
  apf::splitMdsMesh(m, plan, factor, runAfter);
}

void balance(apf::Mesh2* m)
{
  Parma_PrintPtnStats(m, "preRefine");
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  double tolerance = 1.05;
  const double step = 0.2;
  const int verbose = 0;
  apf::Balancer* balancer = Parma_MakeElmBalancer(m, step, verbose);
  balancer->balance(weights, tolerance);
  delete balancer;
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  Parma_PrintPtnStats(m, "postRefine");
}

}

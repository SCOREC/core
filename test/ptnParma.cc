#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <apfZoltan.h>
#include <parma.h>
#include <cassert>
#include <cstdlib>

namespace {

const char* modelFile = 0;
const char* meshFile = 0;
const char* outFile = 0;
const char* method = 0;
const char* approach = 0;
int partitionFactor = 1;
int isLocal = 0;

void setWeight(apf::Mesh* m, apf::MeshTag* tag, int dim) {
  double w = 1;
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

void freeMesh(apf::Mesh* m)
{
  m->destroyNative();
  apf::destroyMesh(m);
}

apf::Migration* getPlan(apf::Mesh* m)
{
  std::map<std::string, int> ztnMethod;
  ztnMethod["rib"] = apf::RIB;
  ztnMethod["rcb"] = apf::RCB;
  ztnMethod["hg"] = apf::HYPERGRAPH;
  ztnMethod["pmetis"] = apf::PARMETIS;

  std::map<std::string, int> ztnApproach;
  ztnApproach["ptn"] = apf::PARTITION;
  ztnApproach["reptn"] = apf::REPARTITION;
  ztnApproach["refine"] = apf::REFINE;
  ztnApproach["kway"] = apf::PART_KWAY;
  ztnApproach["geomkway"] = apf::PART_GEOM_KWAY;
  ztnApproach["refkway"] = apf::REFINE_KWAY;

  apf::Splitter* splitter;
  if( isLocal ) {
    splitter =
      apf::makeZoltanSplitter(m, ztnMethod[method], ztnApproach[approach]);
  } else {
    splitter =
      apf::makeZoltanGlobalSplitter(m, ztnMethod[method], ztnApproach[approach]);
  }
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  apf::Migration* plan = splitter->split(weights, 1.05, partitionFactor);
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  delete splitter;
  return plan;
}

void runAfter(apf::Mesh2* m)
{
  Parma_PrintPtnStats(m, "afterSplit");
  apf::MeshTag* weights = setWeights(m);
  const double step = 0.5; const int verbose = 1;
  apf::Balancer* balancer = Parma_MakeVtxElmBalancer(m, step, verbose);
  balancer->balance(weights, 1.05);
  delete balancer;
  clearTags(m, weights);
  m->destroyTag(weights);
  Parma_PrintPtnStats(m, "final");
  //m->writeNative(outFile);
  freeMesh(m);
}

void getConfig(int argc, char** argv)
{
  if ( argc != 8 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <outMesh> "
             "<factor> <method> <approach> <0:global|1:local>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  meshFile = argv[2];
  outFile = argv[3];
  partitionFactor = atoi(argv[4]);
  method = argv[5];
  approach = argv[6];
  isLocal = atoi(argv[7]);
  if(!PCU_Comm_Self())
    fprintf(stderr, "INPUTS model %s mesh %s out %s factor %d "
       "method %s approach %s\n", modelFile, meshFile, outFile,
       partitionFactor, method, approach);
}

}

int main(int argc, char** argv)
{
  int provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided==MPI_THREAD_MULTIPLE);
  PCU_Comm_Init();
  PCU_Comm_Order(true);
  gmi_register_mesh();
  getConfig(argc,argv);
  apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile);
  Parma_PrintPtnStats(m, "initial");
  splitMdsMesh(m, getPlan(m), partitionFactor, runAfter);
  PCU_Comm_Free();
  MPI_Finalize();
}



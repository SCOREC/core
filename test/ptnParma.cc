#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <lionPrint.h>
#include <apfZoltan.h>
#include <parma.h>
#include <pumi_version.h>
#include <lionPrint.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <cstdlib>
#include <memory>

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

pcu::PCU* getGroupedPCU(pcu::PCU *PCUObj)
{
  int self = PCUObj->Self();
  int groupRank = self / partitionFactor;
  int group = self % partitionFactor;
  MPI_Comm groupComm;
  MPI_Comm_split(MPI_COMM_WORLD, group, groupRank, &groupComm);
  return new pcu::PCU(groupComm);
}

void runParma(apf::Mesh2* m) {
  Parma_PrintPtnStats(m, "afterSplit");
  apf::MeshTag* weights = setWeights(m);
  const double step = 0.5; const int verbose = 1;
  apf::Balancer* balancer = Parma_MakeVtxElmBalancer(m, step, verbose);
  balancer->balance(weights, 1.05);
  delete balancer;
  Parma_PrintPtnStats(m, "postVtxElm");

  // let the percent imbalance increase by 10% e.g. 1.05 -> 1.055
  double elmImb = Parma_GetWeightedEntImbalance(m, weights, m->getDimension());
  double elmImbLimit = 1+((elmImb-1)*1.10);
  balancer = Parma_MakeShapeOptimizer(m, step, verbose);
  balancer->balance(weights, elmImbLimit);
  delete balancer;

  clearTags(m, weights);
  m->destroyTag(weights);
  Parma_PrintPtnStats(m, "final");
}

void mymain(bool ismaster, pcu::PCU *PCUObj)
{
  gmi_model* g = gmi_load(modelFile);
  apf::Mesh2* m = NULL;
  apf::Migration* plan = NULL;
  pcu::PCU *groupedPCUObj = getGroupedPCU(PCUObj);
  if (ismaster) {
    m = apf::loadMdsMesh(modelFile,meshFile,groupedPCUObj);
    Parma_PrintPtnStats(m, "initial");
    plan = getPlan(m);
  }
  //used switchPCU here to load the mesh on the groupedPCU, perform tasks and then call repeatMdsMesh
  //on the globalPCU
  if(m != nullptr) m->switchPCU(PCUObj);
  delete groupedPCUObj;
  m = apf::repeatMdsMesh(m, g, plan, partitionFactor, PCUObj);
  runParma(m);
  m->writeNative(outFile);
  freeMesh(m);
}

void getConfig(int argc, char** argv, pcu::PCU *PCUObj)
{
  if ( argc != 8 ) {
    if ( !PCUObj->Self() )
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
  if(!PCUObj->Self())
    lion_eprint(1, "INPUTS model %s mesh %s out %s factor %d "
       "method %s approach %s isLocal %d\n", modelFile, meshFile, outFile,
       partitionFactor, method, approach, isLocal);
}

}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  auto PCUObj = std::unique_ptr<pcu::PCU>(new pcu::PCU(MPI_COMM_WORLD));
  lion_set_verbosity(1);
  if( !PCUObj.get()->Self() )
    lion_oprint(1, "PUMI version %s Git hash %s\n", pumi_version(), pumi_git_sha());
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  lion_set_verbosity(1);
  getConfig(argc,argv,PCUObj.get());
  if (PCUObj.get()->Self() % partitionFactor)
    mymain(false, PCUObj.get());
  else
    mymain(true, PCUObj.get());
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  }
  MPI_Finalize();
}

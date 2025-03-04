#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <lionPrint.h>
#include <parma.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <apfZoltan.h>
#include <pcu_util.h>
#include <cstdlib>

namespace {

const char* modelFile = 0;
const char* meshFile = 0;
const char* outFile = 0;
int partitionFactor = 1;

void freeMesh(apf::Mesh* m)
{
  m->destroyNative();
  apf::destroyMesh(m);
}

apf::Migration* getPlan(apf::Mesh* m)
{
  apf::Splitter* splitter = apf::makeZoltanSplitter(
      m, apf::GRAPH, apf::PARTITION, false);
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  apf::Migration* plan = splitter->split(weights, 1.05, partitionFactor);
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  delete splitter;
  return plan;
}

void getConfig(int argc, char** argv, pcu::PCU *PCUObj)
{
  if ( argc != 5 ) {
    if ( !PCUObj->Self() )
      printf("Usage: %s <model> <mesh> <outMesh> <factor>\n", argv[0]);
    pcu::Finalize();
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  meshFile = argv[2];
  outFile = argv[3];
  partitionFactor = atoi(argv[4]);
  PCU_ALWAYS_ASSERT(partitionFactor <= PCUObj->Peers());
}

}

int main(int argc, char** argv)
{
  pcu::Init(&argc,&argv);
  {
  pcu::PCU PCUObj;
  lion_set_verbosity(1);
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  getConfig(argc,argv,&PCUObj);
  bool isOriginal = ((PCUObj.Self() % partitionFactor) == 0);
  gmi_model* g = 0;
  g = gmi_load(modelFile);
  apf::Mesh2* m = 0;
  apf::Migration* plan = 0;
  auto groupedPCUObj = PCUObj.Split(
    PCUObj.Self() % partitionFactor, PCUObj.Self() / partitionFactor
  );
  if (isOriginal) {
    m = apf::loadMdsMesh(g, meshFile, groupedPCUObj.get());
    plan = getPlan(m);
  }
  //used switchPCU here to load the mesh on the groupedPCU, perform tasks and then call repeatMdsMesh
  //on the globalPCU
  if(m != nullptr) m->switchPCU(&PCUObj);
  m = apf::repeatMdsMesh(m, g, plan, partitionFactor, &PCUObj);
  Parma_PrintPtnStats(m, "");
  m->writeNative(outFile);
  freeMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  }
  pcu::Finalize();
}

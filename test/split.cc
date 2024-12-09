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
  apf::Splitter* splitter = Parma_MakeRibSplitter(m);
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  apf::Migration* plan = splitter->split(weights, 1.10, partitionFactor);
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

void getConfig(int argc, char** argv, pcu::PCU *PCUObj)
{
  if ( argc != 5 ) {
    if ( !PCUObj->Self() )
      printf("Usage: %s <model> <mesh> <outMesh> <factor>\n", argv[0]);
    MPI_Finalize();
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
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
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
  pcu::PCU *groupedPCUObj = getGroupedPCU(&PCUObj);
  if (isOriginal) {
    m = apf::loadMdsMesh(g, meshFile, groupedPCUObj);
    plan = getPlan(m);
  }
  //used switchPCU here to load the mesh on the groupedPCU, perform tasks and then call repeatMdsMesh
  //on the globalPCU
  if(m != nullptr) m->switchPCU(&PCUObj);
  delete groupedPCUObj;
  m = repeatMdsMesh(m, g, plan, partitionFactor, &PCUObj);
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
  MPI_Finalize();
}


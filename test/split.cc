#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <lionPrint.h>
#include <parma.h>
#ifdef PUMI_HAS_EGADS
#include "gmi_egads.h"
#endif
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

void switchToOriginals()
{
  int self = PCU_Comm_Self();
  int groupRank = self / partitionFactor;
  int group = self % partitionFactor;
  MPI_Comm groupComm;
  MPI_Comm_split(MPI_COMM_WORLD, group, groupRank, &groupComm);
  PCU_Switch_Comm(groupComm);
}

void switchToAll()
{
  MPI_Comm prevComm = PCU_Get_Comm();
  PCU_Switch_Comm(MPI_COMM_WORLD);
  MPI_Comm_free(&prevComm);
  PCU_Barrier();
}

void getConfig(int argc, char** argv)
{
  if ( argc != 5 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <outMesh> <factor>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  meshFile = argv[2];
  outFile = argv[3];
  partitionFactor = atoi(argv[4]);
  PCU_ALWAYS_ASSERT(partitionFactor <= PCU_Comm_Peers());
}

}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
#ifdef PUMI_HAS_EGADS
  gmi_register_egads();
  gmi_egads_start();
#endif
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  getConfig(argc,argv);
  bool isOriginal = ((PCU_Comm_Self() % partitionFactor) == 0);
  gmi_model* g = 0;
  g = gmi_load(modelFile);
  apf::Mesh2* m = 0;
  apf::Migration* plan = 0;
  switchToOriginals();
  if (isOriginal) {
    m = apf::loadMdsMesh(g, meshFile);
    plan = getPlan(m);
  }
  switchToAll();
  m = repeatMdsMesh(m, g, plan, partitionFactor);
  Parma_PrintPtnStats(m, "");
  m->writeNative(outFile);
  freeMesh(m);
#ifdef PUMI_HAS_EGADS
  gmi_egads_stop();
#endif
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  PCU_Comm_Free();
  MPI_Finalize();
}


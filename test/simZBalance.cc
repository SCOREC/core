#include <PCU.h>
#include <ma.h>
#include <MeshSim.h>
#include <SimUtil.h>
#include <SimPartitionedMesh.h>
#include <gmi.h>
#include <gmi_sim.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfSIM.h>
#include <apfZoltan.h>
#include <parma.h>

static void initialize(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  gmi_register_sim();
  SimPartitionedMesh_start(&argc, &argv);
  MS_init();
}

static void finalize() {
  MS_exit();
  SimPartitionedMesh_stop();
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  PCU_Comm_Free();
  MPI_Finalize();
}

static void load_balance(apf::Mesh2* m) {
  Parma_PrintPtnStats(m, "initial");
  ma::Tag* weights = Parma_WeighByMemory(m);
  apf::Balancer* b = makeZoltanBalancer(m, apf::GRAPH, apf::REPARTITION);
  b->balance(weights, 1.10);
  delete b;
  apf::removeTagFromDimension(m, weights, m->getDimension());
  Parma_PrintPtnStats(m, "final");
  m->destroyTag(weights);
}

int main(int argc, char** argv) {
  initialize(argc, argv);
  const char* smd_file = argv[1];
  const char* sms_file = argv[2];
  gmi_model* apf_model = gmi_sim_load(0, smd_file);
  pGModel sim_model = gmi_export_sim(apf_model);
  pParMesh sim_mesh = PM_load(sms_file, sim_model, NULL);
  apf::Mesh2* apf_mesh = apf::createMesh(sim_mesh);
  //apf_mesh->verify(); <- this calls Simmetrix's verify function
  apf::verify(apf_mesh);
  load_balance(apf_mesh);
  apf_mesh->destroyNative();
  apf::destroyMesh(apf_mesh);
  gmi_destroy(apf_model);
  finalize();
}

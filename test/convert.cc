#include <PCU.h>
#include <MeshSim.h>
#include <SimPartitionedMesh.h>
#include <SimUtil.h>
#include <SimParasolidKrnl.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <gmi_sim.h>
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  if (argc != 4) {
    if(0==PCU_Comm_Self())
      std::cerr << "usage: " << argv[0] << " <model file> <simmetrix mesh> <pumi mesh>\n";
    return 0;
  }
  Sim_readLicenseFile(NULL);
  SimPartitionedMesh_start(&argc,&argv);
  SimModel_start();
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  pGModel simModel;
  simModel = GM_load(argv[1], NULL, progress);

  pParMesh sim_mesh = PM_load(argv[2],sthreadNone,0,progress);
  apf::Mesh* simApfMesh = apf::createMesh(sim_mesh);
  
  gmi_register_sim();
  gmi_model* pumiMdl = gmi_import_sim(simModel);
  apf::Mesh2* pumiApfMesh = apf::createMdsMesh(pumiMdl, simApfMesh);
  apf::destroyMesh(simApfMesh);
  M_release(sim_mesh);
  /* partition model ownership update */
  pumiApfMesh->preMigrate_();
  pumiApfMesh->postMigrate_();
  assert(!pumiApfMesh->hasMatching());
  pumiApfMesh->verify();
  pumiApfMesh->writeNative(argv[3]);

  pumiApfMesh->destroyNative();
  apf::destroyMesh(pumiApfMesh);

  Progress_delete(progress);
  SimModel_stop();
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
  PCU_Comm_Free();
  MPI_Finalize();
}

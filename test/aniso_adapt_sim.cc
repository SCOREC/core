#include "aniso_adapt.h"
#include <gmi_sim.h>
#include <apfSIM.h>
#include <SimUtil.h>
#include <SimPartitionedMesh.h>
#include "SimAdvMeshing.h"

int main(int argc, char* argv[])
{
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] <<
      "<model.x_t> <model.smd> <mesh.sms>" << std::endl;
      return 1;
  }
  //Load Mesh
  const char* nativefile = argv[1];
  const char* smdfile = argv[2];
  const char* smsfile = argv[3];
  MPI_Init(&argc, &argv);
  pcu::PCU *PCUObj = new pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  gmi_register_mesh();

  Sim_logOn("anisoadapt.log");
  MS_init();
  SimAdvMeshing_start();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  SimPartitionedMesh_start(&argc,&argv);

  gmi_sim_start();
  gmi_register_sim();

  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);
  gmi_model* mdl_ref = gmi_sim_load(nativefile, smdfile);
  pGModel model = gmi_export_sim(mdl_ref);
  pParMesh mesh = PM_load(smsfile, model, progress);
  ma::Mesh* mesh_ref = apf::createMesh(mesh, PCUObj);

  auto createMeshValues = [mdl_ref, mesh_ref]() 
    { return apf::createMdsMesh(mdl_ref, mesh_ref); };

  adaptTests(createMeshValues);

  M_release(mesh);
  Progress_delete(progress);
  gmi_sim_stop();
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();

  delete PCUObj;
  MPI_Finalize();
}


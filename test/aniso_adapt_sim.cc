#include <iostream>
#include <cstdlib>
#include <filesystem>

#include <lionPrint.h>
#include <pcu_util.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <apfMDS.h>
#include <ma.h>
#include "aniso_adapt.h"

#include <gmi_sim.h>
#include <apfSIM.h>
#include <MeshSim.h>
#include <SimUtil.h>
#include <SimPartitionedMesh.h>
#include <PCU.h>

#include "SimParasolidKrnl.h"
#include "MeshSimAdapt.h"
#include "SimDiscrete.h"
#include "SimAdvMeshing.h"
#include "SimMeshTools.h"

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr << "Usage: " << argv[0] <<
      "<model.x_t> <model.smd> <mesh.sms>" << std::endl;
      return 1;
  }
  //Load Mesh
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
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

  gmi_model* mdl_ref = gmi_sim_load(argv[1], argv[2]);
  pGModel model = gmi_export_sim(mdl_ref);
  pParMesh mesh = PM_load(argv[3], model, progress);
  ma::Mesh* mesh_ref = apf::createMesh(mesh, PCUObj);
  apf::Mesh2* m = apf::createMdsMesh(mdl_ref, mesh_ref);

  //Adapt
  m->verify();
  AnIso sf(m, 3, 1);
  ma::Input* in = ma::makeAdvanced(ma::configure(m, &sf));
  adapt(in);

  //Clean up
  m->destroyNative();
  apf::destroyMesh(m);
  delete PCUObj;
  MPI_Finalize();
}


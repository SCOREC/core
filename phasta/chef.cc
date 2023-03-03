#include <apf.h>
#include <apfMesh.h>
#include <gmi_mesh.h>
#include <PCU.h>
#include <lionPrint.h>
#include <pumi_version.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#include <SimPartitionedMesh.h>
#ifdef HAVE_SIMADVMESHING
#include <SimAdvMeshing.h>
#endif
#endif
#include <pcu_util.h>
#include <chef.h>

/** \file chef.cc
    \brief The Chef command line executable. */

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }
}

/** @brief run the operations requested in "adapt.inp" */
int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_Protect();
  lion_set_verbosity(1);
  if( !PCU_Comm_Self() ) {
    lion_oprint(1,"PUMI Git hash %s\n", pumi_version());
    lion_oprint(1,"PUMI version %s Git hash %s\n", pumi_version(), pumi_git_sha());
  }
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  SimPartitionedMesh_start(0, 0);
#ifdef HAVE_SIMADVMESHING
  SimAdvMeshing_start();
#endif
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  gmi_model* g = 0;
  apf::Mesh2* m = 0;
  std::string inputPath = "adapt.inp";
  if(argc==2) inputPath = argv[1];
  fprintf(stderr, "reading inp %s\n", inputPath.c_str());
  ph::Input in;
  in.load(inputPath.c_str());
  chef::cook(g,m,in);
  freeMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
#ifdef HAVE_SIMADVMESHING
  SimAdvMeshing_stop();
#endif
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  PCU_Comm_Free();
  MPI_Finalize();
}


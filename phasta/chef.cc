#include <apf.h>
#include <apfMesh.h>
#include <gmi_mesh.h>
#include <PCU.h>
#include <pumi_version.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <SimPartitionedMesh.h>
#include <SimAdvMeshing.h>
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
  if( !PCU_Comm_Self() )
    printf("PUMI Git hash %s\n", pumi_version());
#ifdef HAVE_SIMMETRIX
  if( !PCU_Comm_Self() ) printf("Going into from SimUtil_start \n");
  SimUtil_start();
  if( !PCU_Comm_Self() ) printf("Returned from SimUtil_start \n");
  Sim_readLicenseFile(0);
  if( !PCU_Comm_Self() ) printf("Returned from Sim_readLicenseFile \n");
  SimPartitionedMesh_start(0, 0);
  if( !PCU_Comm_Self() ) printf("Returned from SimPartitionedMesh_start \n");
  SimAdvMeshing_start();
  if( !PCU_Comm_Self() ) printf("Returned from SimAdvMeshing_start \n");
  gmi_sim_start();
  if( !PCU_Comm_Self() ) printf("Returned from gmi_sim_start \n");
  gmi_register_sim();
  if( !PCU_Comm_Self() ) printf("Returned from gmi_register_sim \n");
#endif
  gmi_register_mesh();
  gmi_model* g = 0;
  apf::Mesh2* m = 0;
  chef::cook(g,m);
  freeMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  SimPartitionedMesh_stop();
  SimAdvMeshing_stop();
  Sim_unregisterAllKeys();
  SimUtil_stop();
#endif
  PCU_Comm_Free();
  MPI_Finalize();
}


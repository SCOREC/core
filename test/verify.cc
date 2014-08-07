#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <PCU.h>
#include <SimUtil.h>
#include <SimModel.h>

int main(int argc, char** argv)
{
  assert(argc==3);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  Sim_readLicenseFile(0);
  SimModel_start();
  gmi_register_mesh();
  gmi_register_sim();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  m->verify();
  m->destroyNative();
  apf::destroyMesh(m);
  SimModel_stop();
  Sim_unregisterAllKeys();
  PCU_Comm_Free();
  MPI_Finalize();
}




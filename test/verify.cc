#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <PCU.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#endif
#include <cassert>

int main(int argc, char** argv)
{
  assert(argc==3);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
#ifdef HAVE_SIMMETRIX
  SimUtil_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  m->verify();
  m->destroyNative();
  apf::destroyMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimUtil_stop();
#endif
  PCU_Comm_Free();
  MPI_Finalize();
}




#include <apf.h>
#include <apfMesh.h>
#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <PCU.h>
#include <SimUtil.h>
#include <chef.h>

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_Protect();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
  gmi_register_mesh();
  gmi_model* g = 0;
  apf::Mesh2* m = 0;
  chef::cook(g,m);
  freeMesh(m);
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  PCU_Comm_Free();
  MPI_Finalize();
}

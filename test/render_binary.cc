#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <SimUtil.h>
#include <cstdlib>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  SimUtil_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  if ( argc != 4 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <out prefix>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  gmi_register_sim();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  apf::writeBinaryVtkFiles(argv[3], m);
  m->destroyNative();
  apf::destroyMesh(m);
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimUtil_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}


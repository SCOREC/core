#include <apf.h>
#include "apfMDS.h"
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <PCU.h>
#include <SimUtil.h>
#include <SimModel.h>

int main(int argc, char** argv)
{
  assert(argc == 3);
  MPI_Init(&argc,&argv);
  Sim_readLicenseFile(0);
  SimModel_start();
  PCU_Comm_Init();
  PCU_Protect();
  gmi_register_mesh();
  gmi_register_sim();
  //load model and mesh
  double t0 = MPI_Wtime();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  double t1 = MPI_Wtime();
  if (!PCU_Comm_Self())
    std::cout << t1-t0 << " seconds to load\n";
  assert(!alignMdsMatches(m));
  m->verify();
  t0 = MPI_Wtime();
  m->writeNative("out.smb");
  t1 = MPI_Wtime();
  if (!PCU_Comm_Self())
    std::cout << t1-t0 << " seconds to write\n";
  // destroy mds
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  SimModel_stop();
  Sim_unregisterAllKeys();
  MPI_Finalize();
}

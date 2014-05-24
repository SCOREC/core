#include <apf.h>
#include "apfMDS.h"
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <PCU.h>

int main(int argc, char** argv)
{
  assert(argc == 3);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_Protect();
  gmi_register_mesh();
  //load model and mesh
  double t0 = MPI_Wtime();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  double t1 = MPI_Wtime();
  if (!PCU_Comm_Self())
    std::cout << t1-t0 << " seconds to load\n";
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
  MPI_Finalize();
}

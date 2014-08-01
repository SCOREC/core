#include <MeshSim.h>
#include <SimPartitionedMesh.h>
#include <SimUtil.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>

int main(int argc, char** argv)
{
  assert(argc==4);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_sim();
  gmi_register_mesh();

  Sim_readLicenseFile(NULL);
  SimPartitionedMesh_start(&argc,&argv);
  SimModel_start();

  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  apf::writeVtkFiles(argv[3], m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}


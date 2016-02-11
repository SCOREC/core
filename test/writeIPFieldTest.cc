#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <PCU.h>
#include <cassert>

int main(int argc, char** argv)
{
  assert(argc==3);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  m->verify();
  apf::createIPField(m, "Cauchy_Stress", apf::MATRIX, 1);
  apf::writeVtkFiles("mesh_with_IPField", m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

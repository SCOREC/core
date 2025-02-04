#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <lionPrint.h>
#include <pcu_util.h>

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==3);
  MPI_Init(&argc,&argv);
  {
  pcu::PCU pcu_obj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2],&pcu_obj);
  m->verify();
  apf::createIPField(m, "Cauchy_Stress", apf::MATRIX, 1);
  apf::writeVtkFiles("mesh_with_IPField", m);
  m->destroyNative();
  apf::destroyMesh(m);
  }
  MPI_Finalize();
}

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <PCU.h>

int main(int argc, char** argv)
{
  assert(argc == 3);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_null();
  apf::Mesh2* m = apf::loadMdsMesh(".null", argv[1]);
  gmi_model* g = m->getModel();
  gmi_write_dmg(g, argv[2]);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}


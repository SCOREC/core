#include <apfMDS.h>
#include <gmi_null.h>
#include <PCU.h>
#include <apfMesh2.h>
#include <apf.h>
#include <apfMPAS.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_null();
  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 2, false);
  apf::loadMpasMesh(m, argv[1]);
  apf::writeVtkFiles("out", m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

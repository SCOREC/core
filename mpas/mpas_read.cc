#include <apfMDS.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <PCU.h>
#include <apfMesh2.h>
#include <apf.h>
#include <apfMPAS.h>

int main(int argc, char** argv)
{
  if (argc != 4) {
    printf("usage: %s <MPAS NetCDF file> <output model name (.dmg)>"
           " <output mesh name (.smb)>\n", argv[0]);
    return 0;
  }
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_null();
  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 2, false);
  apf::loadMpasMesh(m, argv[1]);
  m->acceptChanges();
  apf::deriveMdsModel(m);
  m->verify();
  gmi_write_dmg(m->getModel(), argv[2]);
  m->writeNative(argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

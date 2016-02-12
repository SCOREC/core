#include <apf.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <cstdlib>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 3 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <in .[b8|lb8].ugrid> <out .smb>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_null();
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m = apf::loadMdsFromUgrid(g, argv[1]);
  apf::deriveMdsModel(m);
  m->verify();
  m->writeNative(argv[2]);
  gmi_write_dmg(g,"mdl.dmg");
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}


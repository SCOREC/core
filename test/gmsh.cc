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
      printf("Usage: %s <in .msh> <out .smb>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_null();
  apf::Mesh2* m = apf::loadMdsFromGmsh(gmi_load(".null"), argv[1]);
  m->verify();
  m->writeNative(argv[2]);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}


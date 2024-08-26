#include <apf.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <lionPrint.h>
#include <cstdlib>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU pcu_obj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  if ( argc != 5 ) {
    if ( !pcu_obj.Self() )
      printf("Usage: %s <in .node> <in .elem> <out .dmg> <out .smb>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_null();
  apf::Mesh2* m = apf::loadMdsFromANSYS(argv[1], argv[2], &pcu_obj);
  m->verify();
  gmi_write_dmg(m->getModel(), argv[3]);
  m->writeNative(argv[4]);
  m->destroyNative();
  apf::destroyMesh(m);
  }
  MPI_Finalize();
}



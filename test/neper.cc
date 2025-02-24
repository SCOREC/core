#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <lionPrint.h>
#include <cstdlib>

int main(int argc, char** argv)
{
  pcu::PCU_Init(&argc,&argv);
  {
  pcu::PCU pcu_obj;
  lion_set_verbosity(1);
  if ( argc != 4 ) {
    if ( !pcu_obj.Self() )
      printf("Usage: %s <model> <in .msh> <out .smb>\n", argv[0]);
    pcu::PCU_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsFromGmsh(gmi_load(argv[1]), argv[2], &pcu_obj);
  m->verify();
  m->writeNative(argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  }
  pcu::PCU_Finalize();
}


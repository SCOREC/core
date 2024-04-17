#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <lionPrint.h>
#include <cstdlib>
#include <memory>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  auto pcu_obj = std::unique_ptr<pcu::PCU>(new pcu::PCU(MPI_COMM_WORLD));
  lion_set_verbosity(1);
  if ( argc != 4 ) {
    if ( !pcu_obj.get()->Self() )
      printf("Usage: %s <model> <in .msh> <out .smb>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsFromGmsh(gmi_load(argv[1]), argv[2], pcu_obj.get());
  m->verify();
  m->writeNative(argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  }
  MPI_Finalize();
}


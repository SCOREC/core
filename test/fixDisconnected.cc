#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <parma.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include <cstdlib>

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc == 4);
  MPI_Init(&argc,&argv);
  {
  pcu::PCU pcu_obj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  if ( argc != 4 ) {
    if ( !pcu_obj.Self() )
      printf("Usage: %s <model> <mesh> <out mesh>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2],&pcu_obj);
  Parma_PrintPtnStats(m, "initial");
  Parma_ProcessDisconnectedParts(m);
  Parma_PrintPtnStats(m, "final");
  m->writeNative(argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  }
  MPI_Finalize();
}

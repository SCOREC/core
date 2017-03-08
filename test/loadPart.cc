#include <apf.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <cstdlib>
#include <pcu_util.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_ALWAYS_ASSERT(PCU_Comm_Peers() == 1);
  if ( argc != 4 ) {
    if ( !PCU_Comm_Self() )
      printf("Load a single part from a partitioned mesh and "
             "write it as a serial part.\n"
             "Usage: %s <in part file> <out mesh file> <out model file (.dmg)>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_null();
  gmi_register_mesh();
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m = apf::loadMdsPart(g, argv[1]);
  apf::deriveMdsModel(m);
  gmi_write_dmg(g, argv[3]);
  m->writeNative(argv[2]);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}




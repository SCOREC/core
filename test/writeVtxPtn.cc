#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <lionPrint.h>
#include <parma.h>
#include <pcu_util.h>
#include <cstdlib>

namespace {
  const char* modelFile = 0;
  const char* meshFile = 0;

  void freeMesh(apf::Mesh* m)
  {
    m->destroyNative();
    apf::destroyMesh(m);
  }

  void getConfig(int argc, char** argv, pcu::PCU *PCUObj)
  {
    if ( argc != 4 ) {
      if ( !PCUObj->Self() )
        printf("Usage: %s <model> <mesh> <out prefix>\n", argv[0]);
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    modelFile = argv[1];
    meshFile = argv[2];
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  gmi_register_mesh();
  getConfig(argc,argv,&PCUObj);
  apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile,&PCUObj);
  Parma_WriteVtxPtn(m,argv[3]);
  freeMesh(m);
  }
  MPI_Finalize();
}

#include <parma.h>
#include <PCU.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <cassert>
#include <cstdlib>

struct GroupCode : public Parma_GroupCode
{
  apf::Mesh2* mesh;
  const char* meshFile;
  void run(int)
  {
    mesh->writeNative(meshFile);
  }
};

int main( int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 5 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <out prefix> <reduction-factor>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  GroupCode code;
  code.mesh = apf::loadMdsMesh(argv[1], argv[2]);
  code.meshFile = argv[3];
  apf::Unmodulo outMap(PCU_Comm_Self(), PCU_Comm_Peers());
  Parma_ShrinkPartition(code.mesh, atoi(argv[4]), code);
  PCU_Comm_Free();
  MPI_Finalize();
}


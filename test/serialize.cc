#include <parma.h>
#include <PCU.h>
#include <apfMDS.h>
#include <gmi_mesh.h>

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
  assert(argc==4);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  GroupCode code;
  code.mesh = apf::loadMdsMesh(argv[1], argv[2]);
  code.meshFile = argv[3];
  apf::Unmodulo outMap(PCU_Comm_Self(), PCU_Comm_Peers());
  Parma_ShrinkPartition(code.mesh, PCU_Comm_Peers(), code);
  PCU_Comm_Free();
  MPI_Finalize();
}


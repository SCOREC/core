#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <parma.h>
#include <cassert>

namespace {
  const char* modelFile = 0;
  const char* meshFile = 0;

  void freeMesh(apf::Mesh* m)
  {
    m->destroyNative();
    apf::destroyMesh(m);
  }

  void getConfig(int argc, char** argv)
  {
    assert(argc==4);
    modelFile = argv[1];
    meshFile = argv[2];
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  getConfig(argc,argv);
  apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile);
  Parma_WriteVtxPtn(m,argv[3]);
  freeMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

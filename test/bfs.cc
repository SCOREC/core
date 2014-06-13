#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <PCU.h>

apf::MeshEntity* grabFirstVertex(apf::Mesh* m)
{
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v = m->iterate(it);
  m->end(it);
  return v;
}

void runBFS(apf::Mesh* m, apf::MeshEntity* startVertex)
{
  /* ... */
}

int main(int argc, char** argv)
{
  assert(argc == 3);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_Protect();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  runBFS(m, grabFirstVertex(m));
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}


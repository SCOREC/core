#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <PCU.h>
#include <malloc.h>

int main(int argc, char** argv)
{
  assert(argc==3);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  fprintf(stderr,"memory before loading\n");
  malloc_stats();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  fprintf(stderr,"memory after loading\n");
  malloc_stats();
  fprintf(stderr,"%d v %d e %d t %d q %d T %d W %d Y\n",
      apf::countEntitiesOfType(m, apf::Mesh::VERTEX),
      apf::countEntitiesOfType(m, apf::Mesh::EDGE),
      apf::countEntitiesOfType(m, apf::Mesh::TRIANGLE),
      apf::countEntitiesOfType(m, apf::Mesh::QUAD),
      apf::countEntitiesOfType(m, apf::Mesh::TET),
      apf::countEntitiesOfType(m, apf::Mesh::PRISM),
      apf::countEntitiesOfType(m, apf::Mesh::PYRAMID));
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}





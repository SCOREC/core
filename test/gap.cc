#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <parma.h>
#include <PCU.h>

int main(int argc, char** argv)
{
  assert(argc == 4);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_Debug_Open();
  gmi_register_mesh();
  //load model and mesh
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  int verbose = 1;
  apf::writeVtkFiles("before",m);
  apf::Balancer* balancer = Parma_MakeShapeOptimizer(m, 0.1, verbose);
  balancer->balance(weights, 1.20);
  delete balancer;
  apf::removeTagFromDimension(m, weights, m->getDimension());
  apf::writeVtkFiles("after",m);

  m->destroyTag(weights);
  m->writeNative(argv[3]);
  // destroy mds
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

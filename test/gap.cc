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
  Parma_PrintPtnStats(m, "initial");
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  int verbose = 2; // set to 1 to silence the 'endStep' stats
  double stepFactor = 0.05;
  apf::Balancer* balancer = Parma_MakeShapeOptimizer(m, stepFactor, verbose);
  balancer->balance(weights, 1.20);
  delete balancer;
  apf::removeTagFromDimension(m, weights, m->getDimension());
  Parma_PrintPtnStats(m, "final");
  m->destroyTag(weights);
  m->writeNative(argv[3]);
  // destroy mds
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

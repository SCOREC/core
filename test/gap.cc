#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <parma.h>
#include <PCU.h>

namespace {
  apf::MeshTag* setWeights(apf::Mesh* m) {
    apf::MeshIterator* it = m->begin(m->getDimension());
    apf::MeshEntity* e;
    apf::MeshTag* tag = m->createDoubleTag("parma_weight", 1);
    double w = 1.0;
    while ((e = m->iterate(it))) 
      m->setDoubleTag(e, tag, &w);
    m->end(it);
    return tag;
  }
}
int main(int argc, char** argv)
{
  assert(argc == 4);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_Comm_Order(true);
  PCU_Debug_Open();
  gmi_register_mesh();
  //load model and mesh
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  apf::writeVtkFiles("before",m);
  Parma_PrintPtnStats(m, "initial");
  apf::MeshTag* weights = setWeights(m);
  int verbose = 2; // set to 1 to silence the 'endStep' stats
  double stepFactor = 0.5;
  apf::Balancer* balancer = Parma_MakeShapeOptimizer(m, stepFactor, verbose);
  balancer->balance(weights, 1.20);
  delete balancer;
  apf::writeVtkFiles("after",m);
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

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <parma.h>
#include <PCU.h>
#include <SimUtil.h>
#include <cassert>

namespace {
  void setWeight(apf::Mesh* m, apf::MeshTag* tag, int dim, double w=1.0) {
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(dim);
    while ((e = m->iterate(it)))
      m->setDoubleTag(e, tag, &w);
    m->end(it);
  }

  apf::MeshTag* setWeights(apf::Mesh* m, double edgeWeight) {
    double weights[4] = {1.,edgeWeight,1.,1.};
    weights[1] = (!PCU_Comm_Self()) ? edgeWeight : 1.0;
    const int d = m->getDimension();
    apf::MeshTag* tag = m->createDoubleTag("parma_weight", 1);
    setWeight(m, tag, 0, weights[0]);
    setWeight(m, tag, 1, weights[1]);
    setWeight(m, tag, d, weights[d]);
    return tag;
  }

  void clearTags(apf::Mesh* m, apf::MeshTag* t) {
    const int d = m->getDimension();
    apf::removeTagFromDimension(m, t, 0);
    apf::removeTagFromDimension(m, t, 1);
    apf::removeTagFromDimension(m, t, d);
  }
}

int main(int argc, char** argv)
{
  assert(argc == 6);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  SimUtil_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  PCU_Comm_Order(true);
  if ( argc != 6 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <out mesh> <edge weight> <tgt imb>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  gmi_register_sim();
  //load model and mesh
  double targetImb = atof(argv[5]);
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  apf::MeshTag* weights = setWeights(m,atof(argv[4]));
  Parma_PrintWeightedPtnStats(m, weights, "initial");
  const double step = 0.5; const int verbose = 1;
  apf::Balancer* balancer = Parma_MakeVtxEdgeElmBalancer(m, step, verbose);
  if( !PCU_Comm_Self() )
    fprintf(stderr, "STATUS target imbalance %.2f\n", targetImb);
  balancer->balance(weights, targetImb);
  delete balancer;
  Parma_PrintWeightedPtnStats(m, weights, "final");
  clearTags(m, weights);
  m->destroyTag(weights);
  m->writeNative(argv[3]);
  // destroy mds
  m->destroyNative();
  apf::destroyMesh(m);
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimUtil_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}

#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <parma.h>
#include <PCU.h>

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
    fprintf(stderr, "STATUS %d entity weights %.2f %.2f %.2f %.2f\n",
        PCU_Comm_Self(), weights[0], weights[1], weights[2], weights[3]);
    apf::MeshTag* tag = m->createDoubleTag("parma_weight", 1);
    for(int i=0; i<=m->getDimension(); i++)
      setWeight(m, tag, i, weights[i]);
    return tag;
  }

  void clearTags(apf::Mesh* m, apf::MeshTag* t) {
    for(int i=0; i<=m->getDimension(); i++)
      apf::removeTagFromDimension(m, t, i);
  }
}

int main(int argc, char** argv)
{
  assert(argc == 5);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_Debug_Open();
  if ( argc != 5 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <out mesh> <edge weight>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  //load model and mesh
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  Parma_PrintPtnStats(m, "initial", true);
  apf::MeshTag* weights = setWeights(m,atof(argv[4]));
  Parma_PrintWeightedPtnStats(m, weights, "initial");
  const double step = 0.5; const int verbose = 1;
  apf::Balancer* balancer = Parma_MakeVtxEdgeElmBalancer(m, step, verbose);
  balancer->balance(weights, 1.10);
  delete balancer;
  Parma_PrintPtnStats(m, "final", true);
  Parma_PrintWeightedPtnStats(m, weights, "final");
  clearTags(m, weights);
  m->destroyTag(weights);
  m->writeNative(argv[3]);
  // destroy mds
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

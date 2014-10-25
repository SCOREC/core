#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <parma.h>
#include <PCU.h>

apf::MeshTag* setVtxWeights(apf::Mesh* m) {
  double w = 1.0;
  apf::MeshTag* tag = m->createDoubleTag("parma_weight", 1);
  //set edge weights
  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) 
    m->setDoubleTag(e, tag, &w);
  m->end(it);
  //set vtx weights
  it = m->begin(0);
  while ((e = m->iterate(it))) 
    m->setDoubleTag(e, tag, &w);
  m->end(it);
  return tag;
}

int main(int argc, char** argv)
{
  assert(argc == 4);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 4 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <out mesh>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  //load model and mesh
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  apf::MeshTag* weights = setVtxWeights(m);
  const double step = 0.1; const int verbose = 1;
  apf::Balancer* balancer = Parma_MakeEdgeBalancer(m, step, verbose);
  balancer->balance(weights, 1.05);
  delete balancer;
  apf::removeTagFromDimension(m, weights, 0);
  apf::removeTagFromDimension(m, weights, 1);
  m->destroyTag(weights);
  m->writeNative(argv[3]);
  // destroy mds
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

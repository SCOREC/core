#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <parma.h>
#include <PCU.h>

apf::MeshTag* setVtxWeights(apf::Mesh* m) {
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;
  apf::MeshTag* tag = m->createDoubleTag("parma_weight", 1);
  double w = 1.0;
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
  gmi_register_mesh();
  //load model and mesh
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  apf::MeshTag* weights = setVtxWeights(m);
  apf::Balancer* balancer = Parma_MakeVtxBalancer(m);
  balancer->balance(weights, 1.05);
  delete balancer;
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  m->writeNative(argv[3]);
  // destroy mds
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

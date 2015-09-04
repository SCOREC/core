#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <parma.h>
#include <PCU.h>
#include <SimUtil.h>
#include <cassert>

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

int main(int argc, char** argv)
{
  assert(argc == 4);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  SimUtil_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  if ( argc != 4 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <out mesh>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  gmi_register_sim();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  double imbalance[4];
  Parma_GetEntImbalance(m,&imbalance);
  if(!PCU_Comm_Self())
    fprintf(stdout, "imbalance <v e f r> %.3f %.3f %.3f %.3f\n",
        imbalance[0], imbalance[1], imbalance[2], imbalance[3]);
  apf::MeshTag* weights = setWeights(m);
  const double step = 0.2; const int verbose = 1;
  apf::Balancer* balancer = Parma_MakeElmBalancer(m, step, verbose);
  balancer->balance(weights, 1.05);
  delete balancer;
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  m->writeNative(argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimUtil_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}

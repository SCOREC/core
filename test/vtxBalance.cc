#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <parma.h>
#include <lionPrint.h>
#ifdef PUMI_HAS_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <pcu_util.h>
#include <cstdlib>

namespace {
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
}

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc == 4);
  pcu::Init(&argc,&argv);
  {
  pcu::PCU PCUObj;
  lion_set_verbosity(1);
  if ( argc != 4 ) {
    if ( !PCUObj.Self() )
      printf("Usage: %s <model> <mesh> <out mesh>\n", argv[0]);
    pcu::Finalize();
    exit(EXIT_FAILURE);
  }
#ifdef PUMI_HAS_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  //load model and mesh
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2],&PCUObj);
  Parma_PrintPtnStats(m, "initial");
  apf::MeshTag* weights = setVtxWeights(m);
  const double step = 0.5; const int verbose = 1;
  apf::Balancer* balancer = Parma_MakeVtxBalancer(m, step, verbose);
  balancer->balance(weights, 1.05);
  delete balancer;
  Parma_PrintPtnStats(m, "final");
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  m->writeNative(argv[3]);
  // destroy mds
  m->destroyNative();
  apf::destroyMesh(m);
#ifdef PUMI_HAS_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  }
  pcu::Finalize();
}

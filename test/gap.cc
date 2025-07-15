#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <gmi_mesh.h>
#include <parma.h>
#include <lionPrint.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <pcu_util.h>
#include <stdlib.h>

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
  PCU_ALWAYS_ASSERT(argc == 5);
  pcu::Init(&argc,&argv);
  {
  pcu::PCU PCUObj;
  lion_set_verbosity(1);
  if ( argc != 5 ) {
    if ( !PCUObj.Self() )
      printf("Usage: %s <model> <mesh> <max elm imb> <out prefix>\n", argv[0]);
    pcu::Finalize();
    exit(EXIT_FAILURE);
  }
#ifdef HAVE_SIMMETRIX
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
  apf::MeshTag* weights = setWeights(m);
  int verbose = 2; // set to 1 to silence the 'endStep' stats
  double stepFactor = 0.5;
  apf::Balancer* balancer = Parma_MakeShapeOptimizer(m, stepFactor, verbose);
  balancer->balance(weights, atof(argv[3]));
  delete balancer;
  apf::removeTagFromDimension(m, weights, m->getDimension());
  Parma_PrintPtnStats(m, "final");
  m->destroyTag(weights);
  m->writeNative(argv[4]);
  // destroy mds
  m->destroyNative();
  apf::destroyMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  }
  pcu::Finalize();
}

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
#include <cstdlib>

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
  PCU_ALWAYS_ASSERT(argc == 4);
#ifndef SCOREC_NO_MPI
  MPI_Init(&argc,&argv);
#else
  (void) argc, (void) argv;
#endif
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  if ( argc != 4 ) {
    if ( !PCUObj.Self() )
      printf("Usage: %s <model> <mesh> <out mesh>\n", argv[0]);
#ifndef SCOREC_NO_MPI
    MPI_Finalize();
#endif
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
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2],&PCUObj);
  double imbalance[4];
  Parma_GetEntImbalance(m,&imbalance);
  if(!m->getPCU()->Self())
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
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  }
#ifndef SCOREC_NO_MPI
  MPI_Finalize();
#endif
}

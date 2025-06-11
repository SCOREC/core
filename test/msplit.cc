#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <lionPrint.h>
#include <parma.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <pcu_util.h>
#include <cstdlib>
#include <apfMETIS.h>

/**
 * \file msplit.cc
 *
 * Test utility for METIS apf::Splitter. Partition mesh from n to n*factor
 * parts.
 */

namespace {

apf::Migration* getPlan(apf::Mesh* m, int factor) {
  apf::Splitter* splitter = apf::makeMETISsplitter(m);
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  double tol = 1.10;
  apf::Migration* plan = splitter->split(weights, tol, factor);
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  delete splitter;
  return plan;
}

struct Args {
  const char *modelFile = nullptr, *meshFile = nullptr, *outFile = nullptr;
  int factor;
};

Args getConfig(int argc, char* argv[], pcu::PCU &PCUObj)
{
  if ( argc != 5 ) {
    if (PCUObj.Self() == 0)
      printf("Usage: %s <model> <mesh> <outMesh> <factor>\n", argv[0]);
    pcu::Finalize();
    exit(EXIT_FAILURE);
  }
  Args args;
  args.modelFile = argv[1];
  args.meshFile = argv[2];
  args.outFile = argv[3];
  args.factor = std::atoi(argv[4]);
  return args;
}

} // namespace

int main(int argc, char* argv[]) {
  lion_set_verbosity(1);
  pcu::Init(&argc,&argv);
  {
  pcu::PCU PCUObj;
  Args args = getConfig(argc, argv, PCUObj);
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  gmi_model* g = 0;
  g = gmi_load(args.modelFile);
  apf::Mesh2* m = 0;
  apf::Migration* plan = nullptr;
  bool isOriginal = (PCUObj.Self() % args.factor == 0);
  auto splitPCU = PCUObj.Split(PCUObj.Self() % args.factor, 0);
  if (isOriginal) {
    m = apf::loadMdsMesh(g, args.meshFile, splitPCU.get());
    plan = getPlan(m, args.factor);
    m->switchPCU(&PCUObj);
  }
  m = repeatMdsMesh(m, g, plan, args.factor, &PCUObj);
  Parma_PrintPtnStats(m, "");
  m->writeNative(args.outFile);
  m->destroyNative();
  apf::destroyMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  MS_exit();
  Sim_unregisterAllKeys();
#endif
  }
  pcu::Finalize();
}


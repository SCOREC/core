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

namespace {

apf::Migration* getPlan(apf::Mesh* m, int factor) {
  apf::Splitter* splitter = apf::makeMETISsplitter(m);
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  apf::Migration* plan = splitter->split(weights, 1.10, factor);
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  delete splitter;
  return plan;
}

struct Args {
  const char *modelFile = nullptr, *meshFile = nullptr, *outFile = nullptr;
  int partitionFactor = 1;
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
  args.partitionFactor = std::atoi(argv[4]);
  PCU_ALWAYS_ASSERT(args.partitionFactor <= PCUObj.Peers());
  return args;
}

} // namespace

int main(int argc, char* argv[]) {
  lion_set_verbosity(1);
  pcu::Init(&argc,&argv);
  {
  pcu::PCU PCUObj;
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  Args args = getConfig(argc, argv, PCUObj);
  bool isOriginal = ((PCUObj.Self() % args.partitionFactor) == 0);
  gmi_model* g = 0;
  g = gmi_load(args.modelFile);
  apf::Mesh2* m = 0;
  apf::Migration* plan = nullptr;
  auto singlePCU = PCUObj.Split(PCUObj.Self(), 0);
  if (isOriginal) {
    m = apf::loadMdsMesh(g, args.meshFile, singlePCU.get());
    plan = getPlan(m, args.partitionFactor);
    m->switchPCU(&PCUObj);
  }
  m = repeatMdsMesh(m, g, plan, args.partitionFactor, &PCUObj);
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


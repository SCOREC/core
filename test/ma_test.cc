#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <gmi.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <apfPartition.h>
#include <lionPrint.h>
#include <parma.h>
#include <apfMETIS.h>
#include <apfZoltan.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <pcu_util.h>
#include <string>
#include <vector>
#include "maLayer.h"

class AnIso : public ma::AnisotropicFunction
{
  public:
    AnIso(ma::Mesh* m, double sf1, double sf2) :
      sizeFactor1(sf1), sizeFactor2(sf2)
    {
      average = ma::getAverageEdgeLength(m);
      ma::getBoundingBox(m, lower, upper);
    }
    virtual void getValue(ma::Entity*, ma::Matrix& R, ma::Vector& H)
    {
      double h = average/sizeFactor1;
      H = ma::Vector(h, h, h/sizeFactor2);
      R = ma::Matrix(
        1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0
      );
    }
  private:
    double sizeFactor1, sizeFactor2, average;
    ma::Vector lower, upper;
};

int main(int argc, char** argv)
{
  bool partitionAtRuntime = false;
  std::vector<const char*> positionalArgs;
  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--partition")
      partitionAtRuntime = true;
    else
      positionalArgs.push_back(argv[i]);
  }
  PCU_ALWAYS_ASSERT(positionalArgs.size() >= 2 && positionalArgs.size() <= 4);
  const char* modelFile = positionalArgs[0];
  const char* meshFile = positionalArgs[1];
  const char* layerTagString =
    (positionalArgs.size() >= 3) ? positionalArgs[2] : "";
  const double adaptRefineFactor =
    (positionalArgs.size() == 4) ? atoi(positionalArgs[3]) : 3;
  pcu::Init(&argc,&argv);
  {
  pcu::PCU PCUObj;
  lion_set_verbosity(1);
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  ma::Mesh* m = nullptr;
  gmi_model* model = nullptr;
  if (partitionAtRuntime && PCUObj.Peers() > 1) {
    model = gmi_load(modelFile);
    auto soloPCU = PCUObj.Split(PCUObj.Self(), 0);
    apf::Migration* plan = nullptr;
    if (PCUObj.Self() == 0) {
      m = apf::loadMdsMesh(model, meshFile, soloPCU.get());
      apf::Splitter* splitter;
#if defined(PUMI_HAS_ZOLTAN)
      splitter = apf::makeZoltanSplitter(m, apf::GRAPH, apf::PARTITION);
#elif defined(PUMI_HAS_METIS)
      splitter = apf::makeMETISsplitter(m);
#else
      splitter = Parma_MakeRibSplitter(m);
#endif
      apf::MeshTag* weights = Parma_WeighByMemory(m);
      plan = splitter->split(weights, 1.10, PCUObj.Peers());
      apf::migrateBoundaryLayerStackTogether(m, plan);
      apf::removeTagFromDimension(m, weights, m->getDimension());
      m->destroyTag(weights);
      delete splitter;
      m->switchPCU(&PCUObj);
    }
    m = apf::repeatMdsMesh(m, model, plan, PCUObj.Peers(), &PCUObj);
    apf::disownMdsModel(m);
  } else {
    m = apf::loadMdsMesh(modelFile, meshFile, &PCUObj);
  }
  m->verify();
  apf::writeVtkFiles("before",m);
  AnIso sf(m, .5, .5);
  ma::Input* in = ma::makeAdvanced(ma::configure(m, &sf));
  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  if(std::string(layerTagString).length()) {
    lion_oprint(1,"disabling adaptation in layer elements tagged %s\n", layerTagString);
    in->userDefinedLayerTagName = layerTagString;
  } else {
    in->shouldSwapLayer = true;
    in->shouldRefineLayer = true;
    in->shouldCoarsenLayer = true;
  }
  ma::adapt(in);
  m->verify();
  m->writeNative("after.smb");
  apf::writeVtkFiles("after",m);
  m->destroyNative();
  apf::destroyMesh(m);
  if (model)
    gmi_destroy(model);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  }
  pcu::Finalize();
}

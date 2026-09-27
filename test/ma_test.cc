#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <lionPrint.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <pcu_util.h>

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
  PCU_ALWAYS_ASSERT(argc>=3);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  const char* layerTagString = (argc==4) ? argv[3] : "";
  const double adaptRefineFactor = (argc==5) ? atoi(argv[4]) : 3;
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
  ma::Mesh* m = apf::loadMdsMesh(modelFile,meshFile,&PCUObj);
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
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  }
  pcu::Finalize();
}


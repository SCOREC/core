#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#include <lionPrint.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <pcu_util.h>

class Linear : public ma::IsotropicFunction
{
  public:
    Linear(ma::Mesh* m)
    {
      mesh = m;
      average = ma::getAverageEdgeLength(m);
      ma::getBoundingBox(m,lower,upper);
    }
    virtual double getValue(ma::Entity* v)
    {
      ma::Vector p = ma::getPosition(mesh,v);
      double x = (p[0] - lower[0])/(upper[0] - lower[0]);
      return average*(4*x+2)/3;
    }
  private:
    ma::Mesh* mesh;
    double average;
    ma::Vector lower;
    ma::Vector upper;
};

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc>=3);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  const char* layerTagString = (argc==4) ? argv[3] : "";
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  ma::Mesh* m = apf::loadMdsMesh(modelFile,meshFile);
  m->verify();
  Linear sf(m);
  ma::Input* in = ma::makeAdvanced(ma::configure(m, &sf));
  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  if(std::string(layerTagString).length()) {
    lion_oprint(1,"disabling adaptation in layer elements tagged %s\n", layerTagString);
    in->userDefinedLayerTagName = layerTagString;
  } else {
    in->shouldRefineLayer = true;
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
  PCU_Comm_Free();
  MPI_Finalize();
}


#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#endif
#include <cassert>

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
  assert(argc==3);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
#ifdef HAVE_SIMMETRIX
  SimUtil_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  ma::Mesh* m = apf::loadMdsMesh(modelFile,meshFile);
  m->verify();
  Linear sf(m);
  ma::Input* in = ma::configure(m, &sf);
  if (!PCU_Comm_Self())
    printf("Matched mesh: disabling"
           " snapping, and shape correction,\n");
  in->shouldSnap = false;
  in->shouldFixShape = false;
  in->shouldTransferParametric = false;
  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->shouldRefineLayer = true;
  ma::adapt(in);
  m->verify();
  apf::writeVtkFiles("after",m);
  m->destroyNative();
  apf::destroyMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimUtil_stop();
#endif
  PCU_Comm_Free();
  MPI_Finalize();
}


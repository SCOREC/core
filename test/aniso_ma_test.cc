#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#include <pcu_util.h>

#include <stdlib.h>

class AnIso : public ma::AnisotropicFunction
{
  public:
    AnIso(ma::Mesh* m)
    {
      mesh = m;
      average = ma::getAverageEdgeLength(m);
      ma::getBoundingBox(m,lower,upper);
    }
    virtual void getValue(ma::Entity* v, ma::Matrix& R, ma::Vector& H)
    {
      ma::Vector p = ma::getPosition(mesh,v);
      double x = (p[0] - lower[0])/(upper[0] - lower[0]);
      double sizeFactor = 2.;
      if (x < 0.5)
	sizeFactor = 3.;
      ma::Vector h(average, average/sizeFactor, average/sizeFactor);
      ma::Matrix r(1.0, 0.0, 0.0,
		   0.0, 1.0, 0.0,
		   0.0, 0.0, 1.0);
      H = h;
      R = r;
    }
  private:
    ma::Mesh* mesh;
    double average;
    ma::Vector lower;
    ma::Vector upper;
};

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==4);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  bool logInterpolation = atoi(argv[3]) > 0 ? true : false;
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  ma::Mesh* m = apf::loadMdsMesh(modelFile,meshFile);
  m->verify();
  apf::writeVtkFiles("aniso_before",m);
  AnIso sf(m);
  ma::Input* in = ma::configure(m, &sf, 0, logInterpolation);
  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->shouldRefineLayer = true;
  in->goodQuality = 0.2;
  ma::adapt(in);
  m->verify();
  if (logInterpolation)
    apf::writeVtkFiles("aniso_log_interpolation_after",m);
  else
    apf::writeVtkFiles("aniso_after",m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}


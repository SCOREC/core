#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <apfField.h>
#include <PCU.h>
#include <lionPrint.h>
#include <pcu_util.h>

#include <stdlib.h>

class AnIso : public ma::AnisotropicFunction
{
  public:
    AnIso(ma::Mesh* m)
    {
      mesh = m;
      ma::getBoundingBox(m,lower,upper);
      targetMetric = m->findField("target_metric");
      PCU_ALWAYS_ASSERT(targetMetric);
      int nComps = targetMetric->countComponents();
      fprintf(stderr, "components %d\n", nComps);
      vals = new double[nComps];
    }
    ~AnIso() {
      delete [] vals;
    }
    virtual void getValue(ma::Entity* v, ma::Matrix& R, ma::Vector& H)
    { 
      apf::getComponents(targetMetric,v,0,vals);
      //convert vals to metric
      ma::Vector h(1.0, 1.0, 1.0);
      ma::Matrix r(1.0, 0.0, 0.0,
		   0.0, 1.0, 0.0,
		   0.0, 0.0, 1.0);
      H = h;
      R = r;
    }
  private:
    ma::Mesh* mesh;
    ma::Vector lower;
    ma::Vector upper;
    apf::Field* targetMetric;
    double* vals;
};

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==4);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  bool logInterpolation = atoi(argv[3]) > 0 ? true : false;
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  gmi_register_mesh();
  gmi_register_null();
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


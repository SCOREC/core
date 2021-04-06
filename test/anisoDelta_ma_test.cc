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
#include "pumi.h"
#include <apfMatrix.h>

#include <stdlib.h>

class AnIso : public ma::AnisotropicFunction {
public:
  AnIso(ma::Mesh* m) {
    mesh = m;
    ma::getBoundingBox(m,lower,upper);
    targetMetric = m->findField("target_metric");
    PCU_ALWAYS_ASSERT(targetMetric);
    int nComps = targetMetric->countComponents();
    vals = new double[nComps];
  }
  ~AnIso() {
    delete [] vals;
  }
  virtual void getValue(ma::Entity* v, ma::Matrix& R, ma::Vector& H) {
    apf::getComponents(targetMetric, v, 0, vals);
    auto scale = 1e4;
    vals[0] = vals[0]/scale;
    vals[1] = vals[1]/scale;
    vals[2] = vals[2]/scale;
    vals[3] = vals[3]/scale;
    vals[4] = vals[4]/scale;
    vals[5] = vals[5]/scale;
    ma::Matrix M(vals[0], vals[3], vals[5],
                 vals[3], vals[1], vals[4],
	         vals[5], vals[4], vals[2]);

    double eigenValues[3];
    ma::Vector eigenVectors[3];
    apf::eigen(M, eigenVectors, eigenValues);
    ma::Matrix RT(eigenVectors); // eigen vectors are stored in the rows of RT
    RT[0] = RT[0].normalize();
    RT[1] = RT[1] - RT[0]*(RT[0]*RT[1]);
    RT[1] = RT[1].normalize();
    RT[2] = apf::cross(RT[0],RT[1]);
    R = apf::transpose(RT);

    double h[3];
    for (int i = 0; i < 3; ++i)
      h[i] = std::sqrt(1.0/(eigenValues[i]*scale));
    H = ma::Vector(h);
  }
private:
  ma::Mesh* mesh;
  ma::Vector lower;
  ma::Vector upper;
  apf::Field* targetMetric;
  double* vals;
};

int main(int argc, char** argv) {
  PCU_ALWAYS_ASSERT(argc == 4);
  const char *modelFile = argv[1];
  const char *meshFile = argv[2];
  bool logInterpolation = false;
  MPI_Init(&argc, &argv);
  PCU_Comm_Init ();
  lion_set_verbosity (1);
  gmi_register_mesh ();
  gmi_register_null ();
  ma::Mesh* m = apf::loadMdsMesh (modelFile, meshFile);
  auto targetMetric = m->findField ("target_metric");
  PCU_ALWAYS_ASSERT (targetMetric);
  m->verify();
  apf::writeVtkFiles ("anisoDelta_before",m);

  AnIso sf (m);
  ma::Input* in = ma::configure (m, &sf, 0, logInterpolation);
  in-> shouldRunPreZoltan = true;
  in-> shouldRunMidParma = true;
  in-> shouldRunPostParma = true;
  in-> shouldRefineLayer = true;
  in-> goodQuality = 0.2;

  ma::adapt (in);
  m-> verify ();

  apf::writeVtkFiles ("anisoDelta_after",m);

  m->destroyNative ();
  apf::destroyMesh (m);
  PCU_Comm_Free ();
  MPI_Finalize ();
}

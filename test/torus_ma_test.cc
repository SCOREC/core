#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#include <lionPrint.h>
#include <pcu_util.h>

class CylindricalShock : public ma::AnisotropicFunction
{
  public:
    CylindricalShock(ma::Mesh* m)
    {
      mesh = m;
      average = ma::getAverageEdgeLength(m);
    }
    virtual void getValue(ma::Entity* v, ma::Matrix& R, ma::Vector& H)
    {
      double cntr = 1.0; // centroid radius for the torus example
      double l = 2.0; // decay factor -- for the torus example
      ma::Vector p = ma::getPosition(mesh,v);
      double x = p[0];
      double y = p[1];
      double norm = std::sqrt(x*x + y*y);

      double distSqr = std::abs(x*x + y*y - cntr*cntr);

      // sizes along principal directions
      ma::Vector h(average * std::abs(1 - exp(-1. * l * distSqr)) / 2. + average / 5.0,
      	           average / 1.,
      	           average / 1.);
      // principal directions
      ma::Matrix r(x/norm,-y/norm, 0.0,
		   y/norm, x/norm, 0.0,
		   0.0,    0.0,    1.0);
      H = h;
      R = r;
    }
  private:
    ma::Mesh* mesh;
    double average;
};

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==3);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  gmi_register_mesh();
  ma::Mesh* m = apf::loadMdsMesh(modelFile,meshFile);
  m->verify();
  apf::writeVtkFiles("torus_before",m);
  CylindricalShock sf(m);
  ma::Input* in = ma::makeAdvanced(ma::configure(m, &sf));
  // Note that the optimal number of iterations will be set
  // inside the call to ma::configure above. However for this
  // test we override this value to 3 to reduce the run time of
  // test.
  in->maximumIterations = 3;
  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  ma::adapt(in);
  m->verify();
  apf::writeVtkFiles("torus_after",m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}


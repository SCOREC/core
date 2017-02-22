#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <PCU.h>
#include <cassert>

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

class CylindricalShock_1 : public ma::AnisotropicFunction
{
  public:
    CylindricalShock_1(ma::Mesh* m)
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
      ma::Vector h(0.4* std::abs(1 - exp(-1. * l * distSqr)) / 2. + 0.4 / 100.0,
      	           0.4/ 1.,
      	           0.4/ 1.);
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

class PlanarShock : public ma::AnisotropicFunction
{
  public:
    PlanarShock(ma::Mesh* m, double inFactor)
    {
      mesh = m;
      average = ma::getAverageEdgeLength(m);
      factor = inFactor;
    }
    virtual void getValue(ma::Entity* v, ma::Matrix& R, ma::Vector& H)
    {
      double loc = 0.75; // location of front
      double l = 1.5; // decay factor for shock
      ma::Vector p = ma::getPosition(mesh,v);
      double x = p[0];

      double distX = std::abs(x - loc);

      /* // sizes along principal directions */
      /* ma::Vector h(average * std::abs(1 - exp(-1. * l * distX)) / 2. + average / 10., */
      /* 	           average / 1., */
      /* 	           average / 1.); */


      // sizes along principal directions
      ma::Vector h(0.4 * std::abs(1 - exp(-1. * l * distX)) + 0.4/factor,
      	           0.4,
      	           0.4);


      /* // principal directions */
      ma::Matrix r(1.0, 0.0, 0.0,
		   0.0, 1.0, 0.0,
		   0.0, 0.0, 1.0);
      H = h;
      R = r;
    }
  private:
    ma::Mesh* mesh;
    double average;
    double factor;
};

class PlanarShock_1 : public ma::AnisotropicFunction
{
  public:
    PlanarShock_1(ma::Mesh* m)
    {
      mesh = m;
      average = ma::getAverageEdgeLength(m);
    }
    virtual void getValue(ma::Entity* v, ma::Matrix& R, ma::Vector& H)
    {
      ma::Vector p = ma::getPosition(mesh,v);
      double y = p[1];
      // sizes along principal directions
      ma::Vector h(0.1,0.001+0.198*std::abs(y-0.5),0.1);
      // principal directions
      ma::Matrix r(1.0, 0.0, 0.0,
		               0.0, 1.0, 0.0,
		               0.0, 0.0, 1.0);
      H = h;
      R = r;
    }
  private:
    ma::Mesh* mesh;
    double average;
};

int main(int argc, char** argv)
{
  assert(argc==3);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  ma::Mesh* m = apf::loadMdsMesh(modelFile,meshFile);
  m->verify();
  apf::writeVtkFiles("torus_before",m);
  CylindricalShock_1 sf(m);
//PlanarShock sf(m,100);
  ma::Input* in = ma::configure(m, &sf);
  in->maximumIterations = 10;
  in->shouldRunPreZoltan = true;
  in->shouldRunMidParma = true;
  in->shouldRunPostParma = true;
  in->shouldRefineLayer = true;
  ma::adapt(in);
  m->verify();
  apf::writeVtkFiles("torus_after",m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}


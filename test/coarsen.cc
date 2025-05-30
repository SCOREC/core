#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfConvert.h>
#include <apf.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include <ma.h>
#include "maAdapt.h"
#include "maCoarsen.h"
#include "apfShape.h"

double vtx_coords[6][3] = {
  {1.0, 0.0, 0.0},
  {1.0, 0.5, 0.0},
  {1.0, 1.0, 0.0},
  {0.5, 0.5, 0.0},
  {0.0, 0.0, 0.0},
  {0.0, 1.0, 0.0}
};

apf::Gid triangles[5][3] = {
  {0, 1, 3},
  {1, 2, 3},
  {0, 3, 4},
  {3, 2, 5},
  {3, 5, 4},
};

class AnIso : public ma::AnisotropicFunction
{
  public:
    AnIso(ma::Mesh* m, double sf1, double sf2) :
      mesh(m), sizeFactor1(sf1), sizeFactor2(sf2)
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
    ma::Mesh* mesh;
    double sizeFactor1, sizeFactor2, average;
    ma::Vector lower, upper;
};


int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  gmi_register_mesh();
  gmi_register_null();

  apf::Gid* conn = &triangles[0][0];
  double* coords = &vtx_coords[0][0];

  int nelem=5;
  int etype=2;
  int nverts=6;
  int dim=2;
 
  //Create mesh
  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, dim, false, &PCUObj);
  apf::GlobalToVert outMap;
  apf::construct(m, conn, nelem, etype, outMap);
  apf::alignMdsRemotes(m);
  apf::deriveMdsModel(m);
  apf::setCoords(m, coords, nverts, outMap);
  outMap.clear();
  m->verify();
  apf::writeVtkFiles("before_coarsen", m);

  //Coarsen
  AnIso sf(m, 0.7, 1);
  ma::Input* in = ma::makeAdvanced(ma::configure(m, &sf));
  in->maximumIterations = 1;
  in->shouldCoarsen = true;
  validateInput(in);
  ma::Adapt* a = new ma::Adapt(in);
  coarsen(a);
  m->verify();

  apf::writeVtkFiles("after_coarsen", m, 1);

  m->destroyNative();
  apf::destroyMesh(m);
  }
  MPI_Finalize();
}


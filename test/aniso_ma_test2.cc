#include <iostream>
#include <cstdlib>

#include <lionPrint.h>
#include <pcu_util.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <apfMDS.h>
#include <ma.h>

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

int main(int argc, char* argv[]) {
  if (argc != 5) {
    std::cerr << "Usage: " << argv[0] <<
      "<model> <mesh> <sizeFactor1> <sizeFactor2>" << std::endl;
      return 1;
  }
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  double sizeFactor1 = std::atof(argv[3]), sizeFactor2 = std::atof(argv[4]);
  MPI_Init(&argc, &argv);
  pcu::PCU *PCUObj = new pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  gmi_register_mesh();
  ma::Mesh* m = apf::loadMdsMesh(modelFile, meshFile, PCUObj);
  m->verify();
  apf::writeVtkFiles("aniso_ma_test2_before",m);
  AnIso sf(m, sizeFactor1, sizeFactor2);
  ma::Input* in = ma::makeAdvanced(ma::configure(m, &sf));
  ma::adapt(in);
  m->verify();
  apf::writeVtkFiles("aniso_ma_test2_after",m);
  m->destroyNative();
  apf::destroyMesh(m);
  delete PCUObj;
  MPI_Finalize();
}


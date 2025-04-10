#include <iostream>
#include <cstdlib>
#include <filesystem>

#include <lionPrint.h>
#include <pcu_util.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <apfMDS.h>
#include <ma.h>
#include "aniso_adapt.h"

int main(int argc, char* argv[]) {
  if (argc != 6) {
    std::cerr << "Usage: " << argv[0] <<
      "<model> <mesh> <sizeFactor1> <sizeFactor2> <outputMesh>" << std::endl;
      return 1;
  }
  //Load Mesh
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  double sizeFactor1 = std::atof(argv[3]), sizeFactor2 = std::atof(argv[4]);
  MPI_Init(&argc, &argv);
  pcu::PCU *PCUObj = new pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  gmi_register_mesh();
  ma::Mesh* m = apf::loadMdsMesh(modelFile, meshFile, PCUObj);

  //Adapt
  refineSnapTest(m);
 
  m->writeNative(argv[5]);

  //Clean up
  m->destroyNative();
  apf::destroyMesh(m);
  delete PCUObj;
  MPI_Finalize();
}


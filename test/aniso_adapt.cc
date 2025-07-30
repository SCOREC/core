#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include <stdlib.h>
#include "aniso_adapt.h"

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==3);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  pcu::Init(&argc,&argv);
  {
  pcu::PCU PCUObj;
  lion_set_verbosity(1);
  gmi_register_mesh();
  ma::Mesh* m = apf::loadMdsMesh(modelFile,meshFile,&PCUObj);
  refineSnapTest(m, 0.5, 1);
  m->destroyNative();
  apf::destroyMesh(m);
  }
  pcu::Finalize();
}


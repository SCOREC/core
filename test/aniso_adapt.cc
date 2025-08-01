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
  coarsenTest(modelFile,meshFile,&PCUObj);
  refineSnapTest(modelFile,meshFile,&PCUObj);
  }
  pcu::Finalize();
}


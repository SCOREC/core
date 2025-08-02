#include "ma.h"
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include <stdlib.h>
#include "aniso_adapt.h"

ma::Mesh* createMesh(const char* modelfile, const char* meshfile, pcu::PCU *PCUObj)
{
  return apf::loadMdsMesh(modelfile,meshfile,PCUObj);
}

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

  auto createMeshValues = [modelFile,meshFile,&PCUObj]() 
    { return createMesh(modelFile,meshFile,&PCUObj); };

  adaptTests(createMeshValues);
  }
  pcu::Finalize();
}


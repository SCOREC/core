#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
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
  ma::Mesh* mesh = apf::loadMdsMesh(modelFile,meshFile,&PCUObj);
  adaptTests(mesh, {0, 0, 0, 5168, 10638, 18363, 25244, 28600, 21121, 6330});
  mesh->destroyNative();
  apf::destroyMesh(mesh);
  }
  pcu::Finalize();
}


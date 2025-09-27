#include <PCU.h>
#include <apfCAP.h>
#include <gmi_cap.h>
#include "aniso_adapt.h"

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);

  if (argc < 2) {
    if(PCUObj.Self()==0)
      std::cerr << "usage: " << argv[0]
        << " <.cre file> or\n"
        << " <.smb file> <model file>\n";
    return EXIT_FAILURE;
  }

  gmi_register_mesh();
  gmi_cap_start();
  gmi_register_cap();
  lion_set_verbosity(1);
  const char* meshFile = argv[1];

  gmi_model* model = gmi_load(meshFile); //Freed in adaptTests
  ma::Mesh* apfCapMesh = apf::createCapMesh(model, &PCUObj);
  apf::disownCapModel(apfCapMesh);
  ma::Mesh* mesh = apf::createMdsMesh(model, apfCapMesh, true);
  apf::disownMdsModel(mesh);
  apf::destroyMesh(apfCapMesh);

  auto createMeshValues = [model, mesh]() 
    { return apf::createMdsMesh(model, mesh, true); };

  adaptTests(createMeshValues);
  apf::destroyMesh(mesh);
  gmi_cap_stop();
  }
  MPI_Finalize();
}
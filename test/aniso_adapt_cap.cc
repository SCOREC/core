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

  gmi_cap_start();
  gmi_register_cap();
  lion_set_verbosity(1);
  const char* meshFile = argv[1];

  gmi_model* model = gmi_load(meshFile);
  ma::Mesh* apfCapMesh = apf::createCapMesh(model, &PCUObj);
  apf::disownCapModel(apfCapMesh);
  ma::Mesh* mesh1 = apf::createMdsMesh(model, apfCapMesh);
  apf::disownMdsModel(mesh1);
  ma::Mesh* mesh2 = apf::createMdsMesh(model, apfCapMesh);
  apf::disownMdsModel(mesh2);

  adaptTests(apfCapMesh, mesh2);

  mesh1->destroyNative();
  mesh2->destroyNative();
  apfCapMesh->destroyNative();
  apf::destroyMesh(mesh1);
  apf::destroyMesh(mesh2);
  apf::destroyMesh(apfCapMesh);
  gmi_cap_stop();
  }
  MPI_Finalize();
}
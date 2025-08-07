#include <PCU.h>
#include <apfCAP.h>
#include <apfMDS.h>
#include <samSz.h>
#include <queue>
#include <ma.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <ma.h>
#include <pcu_util.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include "maAdapt.h"
#include "maCoarsen.h"
#include "maRefine.h"
#include "maSnap.h"
#include "lionPrint.h"
#include "aniso_adapt.h"

ma::Mesh* createMesh(const char* meshFile, pcu::PCU* PCUObj)
{
  gmi_model* model = gmi_load(meshFile);
  ma::Mesh* apfCapMesh = apf::createCapMesh(model, PCUObj);
  apf::disownCapModel(apfCapMesh);
  ma::Mesh* apfMesh = apf::createMdsMesh(model, apfCapMesh, true);
  apf::disownMdsModel(apfMesh);
  apf::destroyMesh(apfCapMesh);
  return apfMesh;
}

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
  gmi_register_null();
  gmi_cap_start();
  gmi_register_cap();
  lion_set_verbosity(1);

  const char* meshFile = argv[1];
  auto createMeshValues = [meshFile, &PCUObj]() 
    { return createMesh(meshFile, &PCUObj); };

  adaptTests(createMeshValues);

  gmi_cap_stop();
  }
  MPI_Finalize();
}
#include <spr.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <PCU.h>

int main(int argc, char** argv)
{
  assert(argc==5);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  const char* outFile = argv[3];
  const int order = atoi(argv[4]);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  apf::Mesh2* mesh = apf::loadMdsMesh(modelFile, meshFile);
  apf::Field* f = 
    apf::createLagrangeField(mesh, "solution", apf::VECTOR, order);
  apf::Field* eps = spr::getGradIPField(f, "eps", order);
  apf::destroyField(f);
  double adaptRatio = 0.1;
  apf::Field* sizef = spr::getSPRSizeField(eps,adaptRatio);
  apf::destroyField(eps);
  writeVtkFiles(outFile,mesh);
  apf::destroyField(sizef);
  mesh->destroyNative();
  apf::destroyMesh(mesh);
  PCU_Comm_Free();
  MPI_Finalize();
}


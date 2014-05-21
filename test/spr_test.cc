#include <apfSPR.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <PCU.h>

int main(int argc, char** argv)
{
  assert(argc==3);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_Protect();
  gmi_register_mesh();
  apf::Mesh2* mesh = apf::loadMdsMesh(modelFile, meshFile);
  apf::changeMeshShape(mesh, apf::getLagrange(2));
  apf::Field* f = apf::createLagrangeField(mesh, "solution", apf::VECTOR, 2);
  apf::Field* eps = apf::getGradIPField(f, "eps", 2);
  apf::destroyField(f);
  double adaptRatio = 0.1;
  apf::Field* sizef = apf::getSPRSizeField(eps,adaptRatio);
  apf::destroyField(eps);
  writeVtkFiles("out",mesh);
  apf::destroyField(sizef);
  mesh->destroyNative();
  apf::destroyMesh(mesh);
  PCU_Comm_Free();
  MPI_Finalize();
}


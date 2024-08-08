#include <spr.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include <cstdlib>
#include <memory>

int main(int argc, char** argv)
{
  const char* exe = argv[0];
  char usage[1024];
  snprintf(usage, 1024, "%s <model.dmg> <in.smb> <out-vtk> <p order>\n", exe);
  PCU_ALWAYS_ASSERT_VERBOSE(argc==5, usage);
  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  const char* outFile = argv[3];
  const int order = atoi(argv[4]);
  MPI_Init(&argc,&argv);
  {
  auto PCUObj = std::unique_ptr<pcu::PCU>(new pcu::PCU(MPI_COMM_WORLD));
  lion_set_verbosity(1);
  gmi_register_mesh();
  apf::Mesh2* mesh = apf::loadMdsMesh(modelFile, meshFile, PCUObj.get());
  if (mesh->findTag("coordinates_edg"))
    mesh->changeShape(apf::getSerendipity(), false);
  apf::Field* f =
    apf::createLagrangeField(mesh, "solution", apf::VECTOR, order);
  apf::Field* eps = spr::getGradIPField(f, "eps", order);
  apf::destroyField(f);
  double adaptRatio = 0.1;
  apf::Field* sizef = spr::getSPRSizeField(eps,adaptRatio,PCUObj.get());
  apf::destroyField(eps);
  writeVtkFiles(outFile,mesh);
  apf::destroyField(sizef);
  mesh->destroyNative();
  apf::destroyMesh(mesh);
  }
  MPI_Finalize();
}


#include <apfMDS.h>
#include <gmi_mesh.h>
#include <apfMPAS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <PCU.h>
#include <sstream>

static std::string makeOutFilename(const char* arg)
{
  if (arg)
    return std::string(arg);
  std::stringstream ss;
  ss << "graph.info.part." << PCU_Comm_Peers();
  return ss.str();
}

int main(int argc, char** argv)
{
  if (argc < 4) {
    printf("usage: %s <model (.dmg)>  <mesh (.smb)> <NetCDF file>"
           " <optional graph filename>\n", argv[0]);
    return 0;
  }
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1], argv[2]);
  std::string outFilename = makeOutFilename(argv[4]);
  apf::writeMpasAssignments(m, argv[3], outFilename.c_str());
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

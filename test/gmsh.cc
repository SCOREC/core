#include <apf.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <lionPrint.h>
#include <cstdlib>
#include <string.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  if ( argc != 5 ) {
    if ( !PCUObj.Self() )
      printf("Usage: %s <in .dmg> <in .msh> <out .smb> <out .dmg>\n"
             "The input .msh and output .smb file names are required. \n"
             "If 'none' is specified as the input model file name then \n"
             "a output model (.dmg) will be written to the specified filename. \n"
             "When a **gmsh v2** .msh is passed in, a minimal model will be created from "
             "the mesh.\n"
             "When a **gmsh v4** .msh is passed in, a topological model will be created "
             "from the geometric model entities defined in the gmsh input file.\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_null();
  gmi_register_mesh();

  std::string model(argv[1]);
  std::string gmsh(argv[2]);
  std::string outMesh(argv[3]);
  std::string outModel(argv[4]);

  const int gmshVersion = apf::gmshMajorVersion(gmsh.c_str());
  fprintf(stderr, "version %d\n", gmshVersion);
  apf::Mesh2* m = NULL;
  if (gmshVersion == 2) {
    if (model.compare("none") == 0) {
      m = apf::loadMdsFromGmsh(gmi_load(".null"), gmsh.c_str(), &PCUObj);
      apf::deriveMdsModel(m);
      gmi_write_dmg(m->getModel(),outModel.c_str());
    } else {
      m = apf::loadMdsFromGmsh(gmi_load(model.c_str()), gmsh.c_str(), &PCUObj);
    }
  } else if (gmshVersion == 4) {
    if (model.compare("none") == 0) {
      m = apf::loadMdsDmgFromGmsh(outModel.c_str(), gmsh.c_str(), &PCUObj);
    }
  }
  m->verify();
  m->writeNative(outMesh.c_str());
  m->destroyNative();
  apf::destroyMesh(m);
  }
  MPI_Finalize();
}


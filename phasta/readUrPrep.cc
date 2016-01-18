#include <apf.h>
#include <apfMesh.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <PCU.h>
#include <chef.h>

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }
  static FILE* openfile_read(ph::Input&, const char* path) {
    return fopen(path, "r");
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_Protect();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  ph::Input ctrl;
  ctrl.load(argv[3]);
  ctrl.openfile_read = openfile_read;
  ctrl.buildMapping = 0; //can't map new vertices from UR
  chef::readAndAttachFields(ctrl,m);
  chef::uniformRefinement(ctrl,m);
  chef::preprocess(m,ctrl);
  freeMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}


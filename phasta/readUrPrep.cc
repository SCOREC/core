#include <apf.h>
#include <apfMesh.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <PCU.h>
#include <chef.h>
#include <memory>

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }
  static FILE* openfile_read(ph::Input&, const char* path, pcu::PCU*) {
    return fopen(path, "r");
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  auto pcu_obj = std::unique_ptr<pcu::PCU>(new pcu::PCU(MPI_COMM_WORLD));
  pcu::Protect();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2], pcu_obj.get());
  ph::Input ctrl;
  ctrl.load(argv[3], pcu_obj.get());
  ctrl.openfile_read = openfile_read;
  ctrl.buildMapping = 0; //can't map new vertices from UR
  chef::readAndAttachFields(ctrl,m);
  chef::uniformRefinement(ctrl,m);
  chef::preprocess(m,ctrl);
  freeMesh(m);
  }
  MPI_Finalize();
}


#include <PCU.h>
#include <chef.h>
#include <gmi_mesh.h>
#include <stdio.h>

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  gmi_register_mesh();
  gmi_model* g = 0;
  apf::Mesh2* m = 0;
  chef::OStream* os = chef::makeOStream();
  chef::cook(g,m,os);
  chef::IStream* is = chef::makeIStream(os);
  chef::destroyOStream(os);
  chef::cook(g,m,is);
  chef::destroyIStream(is);
  freeMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

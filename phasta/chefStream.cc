#include <PCU.h>
#include <chef.h>
#include <phstream.h>
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
  GRStream* grs = makeGRStream();
  chef::cook(g,m,grs);
  RStream* rs = makeRStream();
  attachRStream(grs,rs);
  chef::cook(g,m,rs); //TODO don't re-read the mesh and model
  destroyGRStream(grs);
  destroyRStream(rs);
  freeMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

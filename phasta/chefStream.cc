#include <PCU.h>
#include <chef.h>
#include <phstream.h>
#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <SimUtil.h>
#include <stdio.h>

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }
  void reconfigureChef(ph::Input& ctrl) {
    ctrl.tetrahedronize = 0;
    ctrl.solutionMigration = 1;
    ctrl.outMeshFileName = std::string("bz2:chefStream/");
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
  gmi_register_mesh();
  gmi_model* g = NULL;
  apf::Mesh2* m = NULL;
  GRStream* grs = makeGRStream();
  ph::Input ctrl;
  ctrl.load("adapt.inp");
  chef::cook(g,m,ctrl,grs);
  RStream* rs = makeRStream();
  attachRStream(grs,rs);
  reconfigureChef(ctrl);
  chef::cook(g,m,ctrl,rs);
  destroyGRStream(grs);
  destroyRStream(rs);
  freeMesh(m);
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  PCU_Comm_Free();
  MPI_Finalize();
}

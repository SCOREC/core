/** \file chefStream.cc
    \brief An example use of the chef.h and stream APIs
    \remark This example reads inputs via files, executes
            the operations requested by the "adapt.inp"
            configuration file, and writes PHASTA mesh
            and field data to a stream (using the
            the phstream.h API). The requested Chef operations
            are executed again (the second call to 'cook') using
            the stream inputs (instead of files) and then the
            generates PHASTA files containing the mesh and field
            information.
*/

#include <PCU.h>
#include <chef.h>
#include <phstream.h>
#include <gmi_mesh.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#endif
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
#ifdef HAVE_SIMMETRIX
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
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
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
#endif
  PCU_Comm_Free();
  MPI_Finalize();
}

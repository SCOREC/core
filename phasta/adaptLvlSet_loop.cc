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
#include <chef.h>
#include <phstream.h>
#include <gmi_mesh.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#endif
#include <stdio.h>
#include <unistd.h>

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }
  void reconfigureChef(ph::Input& ctrl) {
    ctrl.restartFileName = "restart";
    ctrl.tetrahedronize = 0;
    ctrl.adaptFlag = 1; //enable adaptation
    ctrl.adaptStrategy = 6; //banded isotropic adaptation around the zero level set
    ctrl.solutionMigration = 1;
    ctrl.outMeshFileName = std::string("bz2:chefStream/");
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  pcu::Protect();
#ifdef HAVE_SIMMETRIX
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  gmi_model* g = NULL;
  apf::Mesh2* m = NULL;
  GRStream* grs = makeGRStream(&PCUObj);
  ph::Input ctrl;
  ctrl.load("adaptLvlSet.inp", &PCUObj);
  //preprocess (define bubbles)
  chef::cook(g,m,ctrl,grs,&PCUObj);
  RStream* rs = makeRStream(&PCUObj);
  attachRStream(grs,rs,&PCUObj);
  reconfigureChef(ctrl);
  // attach the solution field from stream, adapt to the level set,
  //   and then preprocess (redefines bubbles)
  chef::cook(g,m,ctrl,rs,grs,&PCUObj);
  attachRStream(grs,rs,&PCUObj);
  // attach a solution field from stream, adapt to the level set,
  //   and then preprocess (redefines bubbles)
  chef::cook(g,m,ctrl,rs,&PCUObj);
  destroyGRStream(grs,&PCUObj);
  destroyRStream(rs,&PCUObj);
  freeMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
#endif
  }
  MPI_Finalize();
}

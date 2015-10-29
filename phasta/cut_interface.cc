#include "phInterfaceCutter.h"
#include "phAttrib.h"
#include <PCU.h>
#include <SimUtil.h>
#include <gmi_sim.h>
#include <apfMDS.h>
#include <apf.h>
#include <cassert>

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  if (argc != 5) {
    fprintf(stderr,"Usage: %s <model> <attributes> <in mesh> <out mesh>\n", argv[0]);
    return 0;
  }
  PCU_Comm_Init();
  assert(PCU_Comm_Peers() == 1);
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
  char const* modelfile = argv[1];
  char const* attribfile = argv[2];
  char const* meshfile = argv[3];
  char const* outfile = argv[4];
  gmi_model* gm;
  gm = gmi_sim_load(modelfile, attribfile);
  ph::BCs bcs;
  ph::getSimmetrixAttributes(gm, bcs);
  apf::Mesh2* m = apf::loadMdsMesh(gm, meshfile);
  m->verify();
  ph::cutInterface(m, bcs);
  m->verify();
  m->writeNative(outfile);
  m->destroyNative();
  apf::destroyMesh(m);
  gmi_sim_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}

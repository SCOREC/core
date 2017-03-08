#include "phInterfaceCutter.h"
#include "phAttrib.h"
#include <PCU.h>
#include <SimUtil.h>
#include <gmi_sim.h>
#include <apfMDS.h>
#include <apf.h>
#include <pcu_util.h>

char const* modelfile;
char const* attribfile;
char const* meshfile;
char const* outfile;

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  if (argc < 4 || argc > 5) {
    fprintf(stderr,"Usage: %s <model .x_t> <attributes .smd> <in mesh> <out mesh>\n", argv[0]);
    fprintf(stderr,"       to take model and attributes in separate files\n");
    fprintf(stderr,"Usage: %s <model+attributes .smd> <in mesh> <out mesh>\n", argv[0]);
    fprintf(stderr,"       to take combined model and attributes file (by simTranslate)\n");
    return 0;
  }
  PCU_Comm_Init();
  PCU_ALWAYS_ASSERT(PCU_Comm_Peers() == 1);
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
  if (argc == 5) {
    modelfile = argv[1];
    attribfile = argv[2];
    meshfile = argv[3];
    outfile = argv[4];
  }
  else if (argc == 4) {
    attribfile = argv[1];
    meshfile = argv[2];
    outfile = argv[3];
    modelfile = argv[4];
  }
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

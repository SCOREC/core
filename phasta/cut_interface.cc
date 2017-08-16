#include "ph.h"
#include "phInterfaceCutter.h"
#include "phAttrib.h"
#include <PCU.h>
#ifdef HAVE_SIMMETRIX
#include <SimUtil.h>
#include <SimModel.h>
#include <SimPartitionedMesh.h>
#include <gmi_sim.h>
#ifdef HAVE_SIMADVMESHING
#include <SimAdvMeshing.h>
#endif
#endif
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
#ifdef HAVE_SIMMETRIX
  Sim_readLicenseFile(0);
  SimModel_start();
  SimPartitionedMesh_start(0, 0);
#ifdef HAVE_SIMADVMESHING
  SimAdvMeshing_start();
#endif
  gmi_sim_start();
  gmi_register_sim();
#endif
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
  apf::Mesh2* m = ph::loadMesh(gm, meshfile);
  m->verify();
#ifdef HAVE_SIMMETRIX
  if (ph::mesh_has_ext(meshfile, "sms"))
    ph::cutInterfaceSIM(m, bcs);
  else
#endif
    ph::cutInterface(m, bcs);
  m->verify();
  m->writeNative(outfile);
  m->destroyNative();
  apf::destroyMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
#ifdef HAVE_SIMADVMESHING
  SimAdvMeshing_stop();
#endif
  SimPartitionedMesh_stop();
  SimModel_stop();
  Sim_unregisterAllKeys();
#endif
  PCU_Comm_Free();
  MPI_Finalize();
}

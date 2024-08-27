#include "phAttrib.h"
#include "phBC.h"
#include "phInterfaceCutter.h"
#include <ph.h>
#include <lionPrint.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <pcu_util.h>
#include <cstdlib>
#include <iostream>

namespace
{
  const char* modelFile = 0;
  const char* inMesh = 0;
  const char* outMesh = 0;
  char const* attribFile = 0;

void freeMesh(apf::Mesh* m)
{
  m->destroyNative();
  apf::destroyMesh(m);
}

void getConfig(int argc, char** argv, pcu::PCU *pcu_obj)
{
  if (argc < 4 || argc > 5) {
    if ( !pcu_obj->Self() ) {
      lion_eprint(1,"Usage: %s <model .x_t> <attributes .smd> <in mesh> <out mesh>\n", argv[0]);
      lion_eprint(1,"       to take model and attributes in separate files\n");
      lion_eprint(1,"Usage: %s <model+attributes .smd> <in mesh> <out mesh>\n", argv[0]);
      lion_eprint(1,"       to take combined model and attributes file (by simTranslate)\n");}
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  if (argc == 5) {
    modelFile = argv[1];
    attribFile = argv[2];
    inMesh = argv[3];
    outMesh = argv[4];}
  else if (argc == 4) {
    attribFile = argv[1];
    inMesh = argv[2];
    outMesh = argv[3];
    modelFile = argv[4];}
}
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU pcu_obj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  getConfig(argc,argv,&pcu_obj);

  gmi_model* gm;
  gm = gmi_sim_load(modelFile, attribFile);
  apf::Mesh2* m = apf::loadMdsMesh(gm, inMesh, &pcu_obj);
  m->verify();
  ph::BCs bcs;
  ph::getSimmetrixAttributes(gm, bcs);
  if(ph::migrateInterface(m, bcs) == -1)
    ph::fail("no DG interface attributes!");
  m->verify();
  m->writeNative(outMesh);
  apf::writeVtkFiles("test", m);
  freeMesh(m);

#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  }
  MPI_Finalize();
}

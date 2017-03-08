#include "phAttrib.h"
#include "phBC.h"
#include "phInterfaceCutter.h"
#include <ph.h>
#include <PCU.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <SimUtil.h>
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

void getConfig(int argc, char** argv)
{
  if (argc < 4 || argc > 5) {
    if ( !PCU_Comm_Self() ) {
      fprintf(stderr,"Usage: %s <model .x_t> <attributes .smd> <in mesh> <out mesh>\n", argv[0]);
      fprintf(stderr,"       to take model and attributes in separate files\n");
      fprintf(stderr,"Usage: %s <model+attributes .smd> <in mesh> <out mesh>\n", argv[0]);
      fprintf(stderr,"       to take combined model and attributes file (by simTranslate)\n");}
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
  PCU_Comm_Init();
  SimUtil_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_mesh();
  gmi_register_sim();

  getConfig(argc,argv);

  gmi_model* gm;
  gm = gmi_sim_load(modelFile, attribFile);
  apf::Mesh2* m = apf::loadMdsMesh(gm, inMesh);
  m->verify();
  ph::BCs bcs;
  ph::getSimmetrixAttributes(gm, bcs);
  if(ph::migrateInterface(m, bcs) == -1)
    ph::fail("no DG interface attributes!");
  m->verify();
  m->writeNative(outMesh);
  apf::writeVtkFiles("test", m);
  freeMesh(m);

  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimUtil_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}

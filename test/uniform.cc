#include <ma.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <apfMDS.h>
#include <PCU.h>
#include <SimUtil.h>
#include <cassert>
#include <stdlib.h>

const char* modelFile = 0;
const char* meshFile = 0;
const char* outFile = 0;

void getConfig(int argc, char** argv)
{
  if ( argc != 4 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <outMesh>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  meshFile = argv[2];
  outFile = argv[3];
}

int main(int argc, char** argv)
{
  assert(argc==4);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  SimUtil_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_mesh();
  gmi_register_sim();
  getConfig(argc,argv);
  ma::Mesh* m = apf::loadMdsMesh(modelFile,meshFile);
  ma::Input* in = ma::configureUniformRefine(m, 1);
  if (in->shouldSnap) {
    in->shouldSnap = false;
    assert(in->shouldTransferParametric);
  }
  in->shouldFixShape = false;
  ma::adapt(in);
  m->writeNative(outFile);
  m->destroyNative();
  apf::destroyMesh(m);
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimUtil_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}


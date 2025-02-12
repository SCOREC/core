#include <ma.h>
#include <crv.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <lionPrint.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <pcu_util.h>
#include <stdlib.h>

const char* modelFile = 0;
const char* meshFile = 0;
const char* outFile = 0;
int level = 1;

void getConfig(int argc, char** argv, pcu::PCU *PCUObj)
{
  if ( argc != 5 ) {
    if ( !PCUObj->Self() )
      printf("Usage: %s <model> <mesh> <outMesh> <subdivision level>\n", argv[0]);
#ifndef SCOREC_NO_MPI
    MPI_Finalize();
#endif
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  meshFile = argv[2];
  outFile = argv[3];
  level = atoi(argv[4]);
}

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==5);
#ifndef SCOREC_NO_MPI
  MPI_Init(&argc,&argv);
#else
  (void) argc, (void) argv;
#endif
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  getConfig(argc,argv,&PCUObj);
  ma::Mesh* m = apf::loadMdsMesh(modelFile,meshFile,&PCUObj);


  int order = 2;
  crv::BezierCurver bc(m, order, 0);
  bc.run();

  ma::Input* in = ma::makeAdvanced(ma::configureUniformRefine(m, level));
  if (in->shouldSnap) {
    PCU_ALWAYS_ASSERT(in->shouldTransferParametric);
  }
  in->shouldFixShape = false;
  crv::adapt(in);
  crv::writeCurvedVtuFiles(m, apf::Mesh::TRIANGLE, order + 3, outFile);
  crv::writeCurvedWireFrame(m, order + 3, outFile);
  m->destroyNative();
  apf::destroyMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  }
#ifndef SCOREC_NO_MPI
  MPI_Finalize();
#endif
}


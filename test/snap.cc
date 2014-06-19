#include <ma.h>
#include <apf.h>
#include <gmi_sim.h>
#include <apfMDS.h>
#include <PCU.h>
#include <SimUtil.h>
#include <SimModel.h>

int main(int argc, char** argv)
{
  assert(argc==4);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  Sim_readLicenseFile(0);
  SimModel_start();
  gmi_register_sim();
  ma::Mesh* m = apf::loadMdsMesh(argv[1],argv[2]);
  ma::Input* in = ma::configureIdentity(m);
  in->maximumIterations = 1;
  assert(in->shouldSnap);
  ma::adapt(in);
  m->writeNative(argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  SimModel_stop();
  Sim_unregisterAllKeys();
  PCU_Comm_Free();
  MPI_Finalize();
}



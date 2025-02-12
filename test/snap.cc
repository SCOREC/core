#include <ma.h>
#include <apf.h>
#include <gmi_sim.h>
#include <apfMDS.h>
#include <lionPrint.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#include <pcu_util.h>

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==4);
#ifndef SCOREC_NO_MPI
  MPI_Init(&argc,&argv);
#else
  (void) argc, (void) argv;
#endif
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
  ma::Mesh* m = apf::loadMdsMesh(argv[1],argv[2],&PCUObj);
  const ma::Input* in = ma::configureIdentity(m);
  ma::adapt(in);
  m->writeNative(argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
  }
#ifndef SCOREC_NO_MPI
  MPI_Finalize();
#endif
}


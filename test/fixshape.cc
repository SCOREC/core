#include <ma.h>
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

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==4);
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
  ma::Mesh* m = apf::loadMdsMesh(argv[1],argv[2],&pcu_obj);
  ma::Input* in = ma::makeAdvanced(ma::configureIdentity(m));
  in->shouldFixShape = true;
  ma::adapt(in);
  m->writeNative(argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  }
  MPI_Finalize();
}




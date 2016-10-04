#include <crv.h>
#include <crvAdapt.h>

#include <apf.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <PCU.h>
#include <SimUtil.h>
#include <cstdlib>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  SimUtil_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  if ( argc != 5 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <order> <out prefix>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  int order = atoi(argv[3]);
  if(order < 1 || order > 6){
    if ( !PCU_Comm_Self() )
      printf("Only 1st to 6th order supported\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  gmi_register_sim();

  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  crv::BezierCurver bc(m,order,0);
  bc.run();

  ma::Input* in = crv::configureShapeCorrection(m);
  crv::adapt(in);

  m->writeNative(argv[4]);

  m->destroyNative();
  apf::destroyMesh(m);

  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimUtil_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}

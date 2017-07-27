#include <parma.h>
#include <PCU.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <crv.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <pcu_util.h>
#include <cstdlib>

struct GroupCode : public Parma_GroupCode
{
  apf::Mesh2* mesh;
  const char* meshFile;
  void run(int)
  {
    mesh->writeNative(meshFile);
  }
};

int main( int argc, char* argv[])
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 5 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <out prefix> <reduction-factor>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  crv::getBezier(2);//hack to make sure curved meshes can be serialized!
  GroupCode code;
  code.mesh = apf::loadMdsMesh(argv[1], argv[2]);
  code.meshFile = argv[3];
  apf::Unmodulo outMap(PCU_Comm_Self(), PCU_Comm_Peers());
  Parma_ShrinkPartition(code.mesh, atoi(argv[4]), code);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  PCU_Comm_Free();
  MPI_Finalize();
}


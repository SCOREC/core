#include <parma.h>
#include <lionPrint.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
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
  {
  pcu::PCU pcu_obj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  if ( argc != 5 ) {
    if ( !pcu_obj.Self() )
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
  gmi_register_null();
  crv::getBezier(2);//hack to make sure curved meshes can be serialized!
  GroupCode code;
  code.mesh = apf::loadMdsMesh(argv[1], argv[2], &pcu_obj);
  code.meshFile = argv[3];
  apf::Unmodulo outMap(code.mesh->getPCU()->Self(), code.mesh->getPCU()->Peers());
  Parma_ShrinkPartition(code.mesh, atoi(argv[4]), code);
  code.mesh->destroyNative();
  apf::destroyMesh(code.mesh);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  }
  MPI_Finalize();
}


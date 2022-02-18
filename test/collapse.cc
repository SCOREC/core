#include <apf.h>
#include <apfMesh.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <PCU.h>
#include <lionPrint.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#endif
#include <parma.h>
#include <pcu_util.h>
#include <cstdlib>

namespace {

  const char* modelFile = 0;
  const char* meshFile = 0;
  const char* outFile = 0;
  int partitionFactor = 1;

  struct GroupCode : public Parma_GroupCode {
    apf::Mesh2* mesh;
    void run(int) {
      mesh->writeNative(outFile);
    }
  };

  void getConfig(int argc, char** argv)
  {
    if ( argc != 5 ) {
      if ( !PCU_Comm_Self() )
        printf("Usage: mpirun -np <inPartCount> %s <model> <mesh> <outMesh> <factor>\n"
               "Reduce the part count of mesh from inPartCount to inPartCount/factor.\n",
               argv[0]);
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    modelFile = argv[1];
    meshFile = argv[2];
    outFile = argv[3];
    partitionFactor = atoi(argv[4]);
    PCU_ALWAYS_ASSERT(partitionFactor <= PCU_Comm_Peers());
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  PCU_Protect();
#ifdef HAVE_SIMMETRIX
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  getConfig(argc,argv);
  GroupCode code;
  code.mesh = apf::loadMdsMesh(modelFile, meshFile);
  apf::Unmodulo outMap(PCU_Comm_Self(), PCU_Comm_Peers());
  Parma_ShrinkPartition(code.mesh, partitionFactor, code);
  code.mesh->destroyNative();
  apf::destroyMesh(code.mesh);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
#endif
  PCU_Comm_Free();
  MPI_Finalize();
}

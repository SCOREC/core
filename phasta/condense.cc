#include <apf.h>
#include <apfMesh.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <lionPrint.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#endif
#include <chef.h>
#include <parma.h>
#include <pcu_util.h>
#include <cstdlib>

namespace {
  static FILE* openfile_read(ph::Input&, const char* path, pcu::PCU*) {
    return fopen(path, "r");
  }

  struct GroupCode : public Parma_GroupCode {
    apf::Mesh2* mesh;
    ph::Input ctrl;
    void run(int) {
      apf::reorderMdsMesh(mesh,NULL);
      chef::preprocess(mesh,ctrl);
    }
  };

  void checkInputs(int argc, char** argv, pcu::PCU *pcu_obj) {
    if ( argc != 3 ) {
      if ( !pcu_obj->Self() )
        lion_oprint(1,"Usage: %s <control .inp> <reduction-factor>\n", argv[0]);
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if ( !pcu_obj->Self() )
      lion_oprint(1,"Input control file %s reduction factor %s\n", argv[1], argv[2]);
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  {
  pcu::PCU pcu_obj = pcu::PCU(MPI_COMM_WORLD);
  pcu::Protect();
  lion_set_verbosity(1);
  checkInputs(argc,argv,&pcu_obj);
#ifdef HAVE_SIMMETRIX
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  GroupCode code;
  ph::Input in;
  in.load(argv[1], &pcu_obj);
  in.openfile_read = openfile_read;
  code.mesh = apf::loadMdsMesh(in.modelFileName.c_str(), in.meshFileName.c_str(), &pcu_obj);
  chef::readAndAttachFields(in,code.mesh);
  apf::Unmodulo outMap(code.mesh->getPCU()->Self(), code.mesh->getPCU()->Peers());
  code.ctrl = in;
  Parma_ShrinkPartition(code.mesh, atoi(argv[2]), code);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
#endif
  }
  MPI_Finalize();
}

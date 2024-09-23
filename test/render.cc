#include <apf.h>
#include <crv.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <lionPrint.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <cstdlib>
#include <string>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  if ( argc != 4 ) {
    if ( !PCUObj.Self() )
      printf("Usage: %s <model> <mesh> <out prefix>\n", argv[0]);
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

  // This (hack) is here to make sure loadMdsMehs
  // does not fail when the input mesh is curved!
  crv::getBezier(2);

  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2],&PCUObj);

  std::string name = m->getShape()->getName();
  int        order = m->getShape()->getOrder();

  if (name == std::string("Bezier")) {
    crv::writeCurvedVtuFiles(m, apf::Mesh::TRIANGLE, order + 2, argv[3]);
    crv::writeCurvedWireFrame(m, order + 8, argv[3]);
  }
  else
    apf::writeVtkFiles(argv[3], m);
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


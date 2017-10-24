#include <apfMIS.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#endif
#include <cstdlib>
#include <apfShape.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 4 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <out prefix>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
#ifdef HAVE_SIMMETRIX
  SimUtil_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();

  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);

  //We will create a field with the coloring as the values on each element
  apf::Field* coloring = apf::createField(m,"colors",apf::SCALAR,
                                          apf::getConstant(m->getDimension()));

  //Set up a MIS with primary type as elements and adjacencies as vertices.
  apf::MIS* mis = apf::initializeMIS(m,m->getDimension(),0);

  while (apf::getIndependentSet(mis)) {
    //This for loop can be thread parallized safetly
    for (int i=0;i<mis->n;i++) {
      //Independent work can be done here
      apf::setScalar(coloring,mis->ents[i],0,mis->color);
    }
  }
  apf::finalizeMIS(mis);

  //A second coloring over the vertices
  apf::Field* coloring2 = apf::createField(m,"colors2",apf::SCALAR,
                                           m->getShape());
  //Another MIS example where primary types are vertices and
  //    adjacencies are edges
  mis = apf::initializeMIS(m,0,1);

  while (apf::getIndependentSet(mis)) {
    for (int i=0;i<mis->n;i++) {
      apf::setScalar(coloring2,mis->ents[i],0,mis->color);
    }
  }
  apf::finalizeMIS(mis);

  //Vtk files should have a coloring of the elements
  //    to prove the coloring forms independent sets
  apf::writeVtkFiles(argv[3], m);
  m->destroyNative();
  apf::destroyMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimUtil_stop();
#endif
  PCU_Comm_Free();
  MPI_Finalize();
}


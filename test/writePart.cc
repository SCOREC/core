#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <lionPrint.h>
#include <cstdlib>

int main(int argc, char** argv)
{
#ifndef SCOREC_NO_MPI
  MPI_Init(&argc,&argv);
#else
  (void) argc, (void) argv;
#endif
  {
  pcu::PCU pcu_obj;
  lion_set_verbosity(1);
  if ( argc != 5 ) {
    if ( !pcu_obj.Self() )
      printf("Usage: %s <model> <mesh> <part id> <out mesh file>\n", argv[0]);
#ifndef SCOREC_NO_MPI
    MPI_Finalize();
#endif
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1], argv[2], &pcu_obj);
  if (m->getPCU()->Self() == atoi(argv[3]))
    apf::writeMdsPart(m, argv[4]);
  m->destroyNative();
  apf::destroyMesh(m);
  }
#ifndef SCOREC_NO_MPI
  MPI_Finalize();
#endif
}



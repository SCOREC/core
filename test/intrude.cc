#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <lionPrint.h>
#include <maExtrude.h>
#include <cstdlib>
#include <iostream>

int main(int argc, char** argv)
{
#ifndef SCOREC_NO_MPI
  MPI_Init(&argc,&argv);
#else
  (void) argc, (void) argv;
#endif
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  if ( argc != 4 ) {
    if ( !PCUObj.Self() )
      printf("Usage: %s <model> <mesh> <out mesh>\n", argv[0]);
#ifndef SCOREC_NO_MPI
    MPI_Finalize();
#endif
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2],&PCUObj);
  ma::ModelExtrusions extrusions;
  extrusions.push_back(ma::ModelExtrusion(
        m->findModelEntity(1, 2),
        m->findModelEntity(2, 3),
        m->findModelEntity(1, 1)));
  extrusions.push_back(ma::ModelExtrusion(
        m->findModelEntity(2, 2),
        m->findModelEntity(3, 1),
        m->findModelEntity(2, 1)));
  extrusions.push_back(ma::ModelExtrusion(
        m->findModelEntity(0, 2),
        m->findModelEntity(1, 3),
        m->findModelEntity(0, 1)));
  size_t nlayers;
  ma::intrude(m, extrusions, &nlayers);
  if (!m->getPCU()->Self())
    std::cout << "counted " << nlayers << " layers\n";
  m->writeNative(argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  }
#ifndef SCOREC_NO_MPI
  MPI_Finalize();
#endif
}

#include "mpi.h"
#include "pumi.h"

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  pcu::PCU *PCUObj = new pcu::PCU(MPI_COMM_WORLD);
  pGeom g = pumi_geom_load(argv[1], PCUObj, "mesh");
  pMesh m = pumi_mesh_load(g, argv[2], 1, PCUObj);
  pumi_mesh_delete(m);
  pumi_geom_delete(g);
  pumi_finalize(PCUObj);
  MPI_Finalize();
}


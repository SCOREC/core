#include "mpi.h"
#include "pumi.h"

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  pumi_load_pcu(&PCUObj);
  pGeom g = pumi_geom_load(argv[1], "mesh");
  pMesh m = pumi_mesh_load(g, argv[2], 1);
  pumi_mesh_delete(m);
  pumi_geom_delete(g);
  }
  MPI_Finalize();
}


#include "mpi.h"
#include "pumi.h"
#include <memory>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  auto PCUObj = std::unique_ptr<pcu::PCU>(new pcu::PCU(MPI_COMM_WORLD));
  pumi_load_pcu(PCUObj.get());
  pGeom g = pumi_geom_load(argv[1], "mesh");
  pMesh m = pumi_mesh_load(g, argv[2], 1);
  pumi_mesh_delete(m);
  pumi_geom_delete(g);
  }
  MPI_Finalize();
}


#include <iostream>

#include <mpi.h>
#include <PCU.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <apfBox.h>
#include <apfMesh2.h>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  pcu::PCU *PCUObj = new pcu::PCU(MPI_COMM_WORLD);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (rank == 0)
    std::cout << "MPI size: " << size << "; PCU size: " <<
      PCUObj->Peers() << std::endl;
  for (int i = 0; i < size; ++i) {
    if (rank == i)
      std::cout << "Hello from MPI rank: " << rank << "; PCU rank: " <<
        PCUObj->Self() << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
  }
  gmi_register_mesh();
  apf::Mesh2* m = apf::makeMdsBox(1, 1, 1, 1, 1, 1, 0, PCUObj);
  apf::destroyMesh(m);
  delete PCUObj;
  MPI_Finalize();
  return 0;
}

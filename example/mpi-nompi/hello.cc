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
  // Print rank info.
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
  // Test SCOREC functions.
  gmi_register_mesh();
  apf::Mesh2* m = apf::makeMdsBox(1, 1, 1, 1, 1, 1, 0, PCUObj);
  apf::destroyMesh(m);
  // MPI Allreduce
  int val = rank + 1, sum;
  MPI_Allreduce(&val, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  int pval = rank + 1, psum;
  psum = PCUObj->Add(pval);
  for (int i = 0; i < size; ++i) {
    if (rank == i)
      std::cout << "MPI(" << rank << ") sum: " << sum <<
      "; PCU sum: " << psum << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
  }
  delete PCUObj;
  MPI_Finalize();
  return 0;
}

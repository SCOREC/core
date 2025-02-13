#include <iostream>

#include <mpi.h>
#include <PCU.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <apfBox.h>
#include <apfMesh2.h>

void print_array(int *arr, size_t n);

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
  // Test Allgather
  int *vals1 = new int[size];
  MPI_Allgather(&val, 1, MPI_INT, vals1, 1, MPI_INT, MPI_COMM_WORLD);
  int *pvals1 = new int[PCUObj->Peers()];
  PCUObj->Allgather(&val, pvals1, 1);
  for (int i = 0; i < size; ++i) {
    if (rank == i) {
      std::cout << "MPI(" << rank << ") gather: ";
      print_array(vals1, size);
      std::cout << "; PCU gather: ";
      print_array(pvals1, PCUObj->Peers());
      std::cout << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  delete vals1;
  delete pvals1;
  int val2[2] = {val, -val}, *vals2 = new int[2 * size];
  MPI_Allgather(val2, 2, MPI_INT, vals2, 2, MPI_INT, MPI_COMM_WORLD);
  int *pvals2 = new int[2 * PCUObj->Peers()];
  PCUObj->Allgather(val2, pvals2, 2);
  for (int i = 0; i < size; ++i) {
    if (rank == i) {
      std::cout << "MPI(" << rank << ") gather2: ";
      print_array(vals2, 2 * size);
      std::cout << "; PCU gather2: ";
      print_array(pvals2, 2 * PCUObj->Peers());
      std::cout << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  delete vals2;
  delete pvals2;
  delete PCUObj;
  MPI_Finalize();
  return 0;
}

void print_array(int *arr, size_t n) {
  std::cout << "[";
  for (int i = 0; i < n; ++i) {
    if (i != 0) std::cout << ", ";
    std::cout << arr[i];
  }
  std::cout << "]";
}

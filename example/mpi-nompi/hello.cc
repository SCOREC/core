#include <iostream>
#include <vector>

#include <mpi.h>
#include <PCU.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <apfBox.h>
#include <apfMesh2.h>

template <typename T>
void print_vector(const std::vector<T> &v);

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  { // pcu object scope
  pcu::PCU PCUObj;
  // Print rank info.
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (rank == 0)
    std::cout << "MPI size: " << size << "; PCU size: " <<
      PCUObj.Peers() << std::endl;
  for (int i = 0; i < size; ++i) {
    if (rank == i)
      std::cout << "Hello from MPI rank: " << rank << "; PCU rank: " <<
        PCUObj.Self() << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // Test SCOREC functions.
  gmi_register_mesh();
  apf::Mesh2* m = apf::makeMdsBox(1, 1, 1, 1, 1, 1, 0, &PCUObj);
  apf::destroyMesh(m);
  // MPI Allreduce
  int val = rank + 1, sum;
  MPI_Allreduce(&val, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  int pval = rank + 1, psum;
  psum = PCUObj.Add(pval);
  for (int i = 0; i < size; ++i) {
    if (rank == i)
      std::cout << "MPI(" << rank << ") sum: " << sum <<
      "; PCU sum: " << psum << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
  }
  // Test Allgather
  std::vector<int> vals1(size);
  MPI_Allgather(&val, 1, MPI_INT, vals1.data(), 1, MPI_INT, MPI_COMM_WORLD);
  std::vector<int> pvals1(PCUObj.Peers());
  PCUObj.Allgather(&val, pvals1.data(), 1);
  for (int i = 0; i < size; ++i) {
    if (rank == i) {
      std::cout << "MPI(" << rank << ") gather: ";
      print_vector(vals1);
      std::cout << "; PCU gather: ";
      print_vector(pvals1);
      std::cout << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  int val2[2] = {val, -val};
  std::vector<int> vals2(2 * size);
  MPI_Allgather(val2, 2, MPI_INT, vals2.data(), 2, MPI_INT, MPI_COMM_WORLD);
  std::vector<int> pvals2(2 * PCUObj.Peers());
  PCUObj.Allgather(val2, pvals2.data(), 2);
  for (int i = 0; i < size; ++i) {
    if (rank == i) {
      std::cout << "MPI(" << rank << ") gather2: ";
      print_vector(vals2);
      std::cout << "; PCU gather2: ";
      print_vector(pvals2);
      std::cout << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  } // pcu object scope
  MPI_Finalize();
  return 0;
}

template <typename T>
void print_vector(const std::vector<T>& v) {
  std::cout << "[";
  for (int i = 0; i < v.size(); ++i) {
    if (i != 0) std::cout << ", ";
    std::cout << v[i];
  }
  std::cout << "]";
}

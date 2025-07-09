#include <iostream>
#include <vector>

#include <mpi.h>
#include <PCU.h>
#include <pcu_util.h>

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
  #ifdef PUMI_NO_MPI
  PCU_ALWAYS_ASSERT(PCUObj.Self() == 0);
  PCU_ALWAYS_ASSERT(PCUObj.Peers() == 1);
  #else
  PCU_ALWAYS_ASSERT(PCUObj.Self() == rank);
  PCU_ALWAYS_ASSERT(PCUObj.Peers() == size);
  #endif
  // MPI Allreduce
  int val = rank + 1, sum;
  MPI_Allreduce(&val, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  int pval = PCUObj.Self() + 1, psum;
  psum = PCUObj.Add(pval);
  for (int i = 0; i < size; ++i) {
    PCU_ALWAYS_ASSERT(sum == (size * (size + 1)) / 2);
    PCU_ALWAYS_ASSERT(psum == (PCUObj.Peers() * (PCUObj.Peers() + 1)) / 2);
  }
  // Test Allgather
  std::vector<int> vals1(size);
  MPI_Allgather(&val, 1, MPI_INT, vals1.data(), 1, MPI_INT, MPI_COMM_WORLD);
  std::vector<int> pvals1(PCUObj.Peers());
  PCUObj.Allgather(&pval, pvals1.data(), 1);
  for (int i = 0; i < size; ++i)
    PCU_ALWAYS_ASSERT(vals1[i] == i + 1);
  for (int i = 0; i < PCUObj.Peers(); ++i)
    PCU_ALWAYS_ASSERT(pvals1[i] == i + 1);
  int val2[2] = {val, -val};
  std::vector<int> vals2(2 * size);
  MPI_Allgather(val2, 2, MPI_INT, vals2.data(), 2, MPI_INT, MPI_COMM_WORLD);
  int pval2[2] = {pval, -pval};
  std::vector<int> pvals2(2 * PCUObj.Peers());
  PCUObj.Allgather(pval2, pvals2.data(), 2);
  for (int i = 0; i < size; ++i) {
    PCU_ALWAYS_ASSERT(vals2[i * 2] == i + 1);
    PCU_ALWAYS_ASSERT(vals2[i * 2 + 1] == -(i + 1));
  }
  for (int i = 0; i < PCUObj.Peers(); ++i) {
    PCU_ALWAYS_ASSERT(pvals2[i * 2] == i + 1);
    PCU_ALWAYS_ASSERT(pvals2[i * 2 + 1] == -(i + 1));
  }
  } // pcu object scope
  MPI_Finalize();
  return 0;
}


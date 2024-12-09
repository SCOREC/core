#include <PCU.h>
#include <pcu_io.h> //pcu_swap_doubles
#include <pcu_util.h> //PCU_ALWAYS_ASSERT
#include <numeric> //iota
#include <iostream> //cerr

int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  const size_t n = 2;
  double *d_orig = new double[n];
  std::iota(d_orig,d_orig+n,0);
  double *d = new double[n];
  std::iota(d,d+n,0);
  pcu_swap_doubles(d, n);
  pcu_swap_doubles(d, n);
  for(size_t i=0; i<n; i++) {
    if(d[i] != d_orig[i]) {
      std::cerr << "ERROR: swap does not match"
                << " d[" << i << "] = " << d[i]
                << " d_orig[" << i << "] = " << d_orig[i] << "\n";

    }
    PCU_ALWAYS_ASSERT(d[i]==d_orig[i]);
  }
  delete [] d_orig;
  delete [] d;
  }
  MPI_Finalize();
  return 0;
}

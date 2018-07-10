#include <PCU.h>
#include <pcu_util.h>

int main(int ac, char * av[])
{
  MPI_Init(&ac,&av);
  int mpi_sz = 0;
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_sz);
  int mpi_rnk = -1;
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_rnk);
  PCU_Comm_Init();
  int sz = PCU_Comm_Peers();
  PCU_ALWAYS_ASSERT(sz == mpi_sz);
  int rnk = PCU_Comm_Self();
  PCU_ALWAYS_ASSERT(rnk == mpi_rnk);
  int * gather_bfr = new int[sz];
  int * compar_bfr = new int[sz];
  for(int rr = 0; rr < sz; ++rr)
    compar_bfr[rr] = rr;
  PCU_Allgather_Int(rnk,gather_bfr);
  for(int rr = 0; rr < sz; ++rr)
    PCU_ALWAYS_ASSERT(gather_bfr[rr] == compar_bfr[rr]);
  delete [] compar_bfr;
  delete [] gather_bfr;
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}

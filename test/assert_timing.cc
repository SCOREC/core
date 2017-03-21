#include <PCU.h>
#include <pcu_util.h>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>

double check_c_assert() {
  double t0 = PCU_Time();
  for (int i = 0; i < 100000000; ++i) {
    double j = (double)i;
    double k = (double)i + 1.0;
    assert(pow(j, 1.01) < pow(k, 1.02));
  }
  double t1 = PCU_Time();
  return t1-t0;
}

double check_pcu_assert() {
  double t0 = PCU_Time();
  for (int i = 0; i < 100000000; ++i) {
    double j = (double)i;
    double k = (double)i + 1.0;
    PCU_ALWAYS_ASSERT(pow(j, 1.01) < pow(k, 1.02));
  }
  double t1 = PCU_Time();
  return t1-t0;
}

int main(int argc, char** argv) {
  PCU_ALWAYS_ASSERT(argc == 2);
  int opt = atoi(argv[1]);
  MPI_Init(0,0);
  PCU_Comm_Init();
  /* i'm avoiding conditionals inside for loops b/c
     i'm paranoid about the timings even though timings
     should not be affected by them at all... */
  if (opt == 0)
    for (int i = 0; i < 5; ++i)
      printf("c assert in %f seconds\n", check_c_assert());
  else
    for (int i = 0; i < 5; ++i)
      printf("pcu assert in %f seconds\n", check_pcu_assert());
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}

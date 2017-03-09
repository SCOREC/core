#include <PCU.h>
#include <pcu_util.h>
#include <cstdio>
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
  printf("pcu assert in %f seconds\n", t1-t0);
}

int main() {
  MPI_Init(0,0);
  PCU_Comm_Init();
  for (int i = 0; i < 5; ++i) {
    double c_time = check_c_assert();
    printf("c assert in %f seconds\n", c_time);
    double pcu_time = check_pcu_assert();
    printf("pcu assert in %f seconds\n", pcu_time);
  }
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}

#include <PCU.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>

double check_c_assert() {
  double t0 = pcu::Time();
  for (int i = 0; i < 100000000; ++i) {
    assert(pow((double)i, 1.01) < pow((double)i +1.0, 1.02));
  }
  double t1 = pcu::Time();
  return t1-t0;
}

double check_pcu_assert() {
  double t0 = pcu::Time();
  for (int i = 0; i < 100000000; ++i) {
    double j = (double)i;
    double k = (double)i + 1.0;
    PCU_ALWAYS_ASSERT(pow(j, 1.01) < pow(k, 1.02));
  }
  double t1 = pcu::Time();
  return t1-t0;
}

int main(int argc, char** argv) {
  PCU_ALWAYS_ASSERT(argc == 2);
  int opt = atoi(argv[1]);
  MPI_Init(0,0);
  {
  pcu::PCU pcu_obj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  /* i'm avoiding conditionals inside for loops b/c
     i'm paranoid about the timings even though timings
     should not be affected by them at all... */
  if (opt == 0)
    for (int i = 0; i < 5; ++i)
      printf("c assert in %f seconds\n", check_c_assert());
  else
    for (int i = 0; i < 5; ++i)
      printf("pcu assert in %f seconds\n", check_pcu_assert());
  }
  MPI_Finalize();
  return 0;
}

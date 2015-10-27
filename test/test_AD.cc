#include "mthAD.h"
#include <cassert>
#include <cstdio>

typedef mth::AD<3> ScalarT;

void print(ScalarT x)
{
  printf("%f [ ", x.val());
  for (unsigned int i=0; i < x.size(); ++i)
    printf(" %f ", x.dx(i));
  printf("]\n");
}

void compare(ScalarT x, double val, double da, double db, double dc)
{
  assert(fabs(x.val() - val) < 1e-15);
  assert(fabs(x.dx(0) - da) < 1e-15);
  assert(fabs(x.dx(1) - db) < 1e-15);
  assert(fabs(x.dx(2) - dc) < 1e-15);
}

/*void compare(ScalarT x, const char* name, const char* val)
{
  printf("%s should be : %s\n", name, val);
  printf("%s is        : ", name);
  print(x);
}*/

ScalarT f1(ScalarT a, ScalarT b, ScalarT c)
{
  return a + b + c;
}

ScalarT f2(ScalarT a, ScalarT b, ScalarT c)
{
  return -a - b - c;
}

ScalarT f3(ScalarT b, ScalarT c)
{
  return b * exp(c);
}

ScalarT f4(ScalarT a, ScalarT c)
{
  return a/c;
}

int main()
{
  ScalarT a = 1.0;
  ScalarT b = 2.0;
  ScalarT c = 3.0;

  a.diff(0);
  b.diff(1);
  c.diff(2);

  ScalarT f = f1(a,b,c);
  compare(f, 6.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000);

  ScalarT g = f2(a,b,c);
  compare(g, -6.000000000000000, -1.000000000000000, -1.000000000000000, -1.000000000000000);

  ScalarT h = f3(b,c);
  compare(h, 40.171073846375336, 0.000000000000000, 20.085536923187668, 40.171073846375336);

  ScalarT i = f4(a,c);
  compare(i, 0.333333333333333, 0.333333333333333, 0.000000000000000, -0.111111111111111);

  ScalarT j = 2.0 * (a + b - c) / 4.0;
  compare(j, 0.000000000000000, 0.500000000000000, 0.500000000000000, -0.500000000000000);

  ScalarT k = (f + h) * i * exp(j*g) / i;
  compare(k, 46.171073846375336, -137.513221539126022, -117.427684615938347, 179.684295385501315);
}


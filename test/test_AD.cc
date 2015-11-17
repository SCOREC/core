#include "mthAD.h"
#include <cassert>

typedef mth::AD<double, 3> AD;

void compare(AD const& x, double val, double da, double db)
{
  assert(fabs(x.val() - val) < 1e-15);
  assert(fabs(x.dx(0) - da) < 1e-15);
  assert(fabs(x.dx(1) - db) < 1e-15);
}

AD f1(AD const& a, AD const& b)
{
  return a + b;
}

int main()
{
  AD a = 1.0;
  AD b = 2.0;

  a.diff(0);
  b.diff(1);

  AD f = f1(a,b);
  compare(f, 3.000000000000000, 1.000000000000000, 1.000000000000000);
}

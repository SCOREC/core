#include "mthQR.h"
#include <iostream>

int main()
{
  mth::Matrix3x3<double> a(
      1, 1, 0,
      0, 1, 0,
      0, 0, 0);
  std::cout << "A\n" << a;
  mth::Matrix<double,3,3> q;
  mth::Matrix<double,3,3> r;
  unsigned rank = decomposeQR(a, q, r);
  std::cout << "rank " << rank << '\n';
  std::cout << "Q\n" << q;
  std::cout << "R\n" << r;
  std::cout << "Q*R\n" << q * r;
}

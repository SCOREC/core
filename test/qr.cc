#include "mthQR.h"
#include "mth_def.h"
#include <iostream>

int main()
{
  double a_dat[3][3] = {
    {1,1,0},
    {0,1,0},
    {0,0,0}
  };
  mth::Matrix<double,0,0> a(3,3);
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j)
    a(i,j) = a_dat[i][j];
  std::cout << "A\n" << a;
  mth::Matrix<double,0,0> q;
  mth::Matrix<double,0,0> r;
  unsigned rank = decomposeQR(a, q, r);
  std::cout << "rank " << rank << '\n';
  std::cout << "Q\n" << q;
  std::cout << "R\n" << r;
  mth::Matrix<double,0,0> qr;
  mth::multiply(q, r, qr);
  std::cout << "Q*R\n" << qr;
}

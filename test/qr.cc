#include "mthQR.h"
#include "mth_def.h"
#include <iostream>

int main()
{
  double a_dat[3][3] = {
    {1,0,0},
    {0,1,0},
    {0,0,1}
  };
  double b_dat[3] = {
    5,
    6,
    7
  };
  mth::Matrix<double,0,0> a(3,3);
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j)
    a(i,j) = a_dat[i][j];
  std::cout << "A\n" << a;
  mth::Vector<double,0> b(3);
  for (unsigned i = 0; i < 3; ++i)
    b(i) = b_dat[i];
  std::cout << "B\n" << b;
  mth::Vector<double,0> x;
  mth::solveQR(a, b, x);
  std::cout << "X\n" << x;
}

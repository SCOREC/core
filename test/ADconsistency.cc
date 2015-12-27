#include <mthAD.h>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <string>
using namespace mth;
void printAndCompare(AD<double, 0> dyn, AD<double, 3> sta, std::string type)
{
  std::cout << "Testing output values for " << type << '\n';
  std::cout << "Static value is: \t\t" << sta.val() << '\n';
  std::cout << "Dynamic value is: \t\t" << dyn.val() << '\n';
  assert(dyn.val() - sta.val() < 1e-14);
  for(unsigned int i = 0; i < 3; i++)
  {
    std::cout << i << " static derivative is: \t" << sta.dx(i) << '\n';
    std::cout << i << " dynamic derivative is:\t" << dyn.dx(i) << '\n';
    assert(dyn.dx(i) - sta.dx(i) < 1e-14);
  }
}
template <class T, unsigned int N>
AD<T, N> addition(AD<T, N> a, AD<T, N> b, AD<T, N> c)
{
  return a + b + c + 37;
}

template <class T, unsigned int N>
AD<T, N> subtraction(AD<T, N> a, AD<T, N> b, AD<T, N> c)
{
  return a - b - c - 37;
}

template <class T, unsigned int N>
AD<T, N> multiplication(AD<T, N> a, AD<T, N> b, AD<T, N> c)
{
  return a * b * c * 37;
}

template <class T, unsigned int N>
AD<T, N> division(AD<T, N> a, AD<T, N> b, AD<T, N> c)
{
  return a / b / c;
}

template <class T, unsigned int N>
AD<T, N> powTest(AD<T, N> a, AD<T, N> b)
{
  return pow(a, b);
}

template <class T, unsigned int N>
AD<T, N> powTest2(AD<T, N> a)
{
  return pow(a, 3.0);
}

template <class T, unsigned int N>
AD<T, N> powTest3(AD<T, N> a)
{
  return pow(3.0, a);
}

template <class T, unsigned int N>
AD<T, N> powTest4(AD<T, N> a)
{
  return pow(3, sin(a));
}

template <class T, unsigned int N>
AD<T, N> complex(AD<T, N> a, AD<T, N> b, AD<T, N> c)
{
  return pow(a, 12) * sqrt(b / c) / sin(c * b) * 37 / c;
}

int main()
{
  std::ios::fmtflags f( std::cout.flags() );
  std::cout << std::setprecision(20);
  AD<double, 3> a = 1.0;
  AD<double, 3> b = 3.0;
  AD<double, 3> c = 7.0;
  a.diff(0);
  b.diff(1);
  c.diff(2);
  AD<double, 0> x = 1.0;
  AD<double, 0> y = 3.0;
  AD<double, 0> z = 7.0;
  x.diff(0, 3);
  y.diff(1, 3);
  z.diff(2, 3);
  printAndCompare(x, a, "sanity check");
  printAndCompare(y, b, "sanity check");
  printAndCompare(z, c, "sanity check");
  AD<double, 3> sta = addition(a, b, c);
  AD<double, 0> dyn = addition(x, y, z);
  printAndCompare(dyn, sta, "addition");
  sta = subtraction(a, b, c);
  dyn = subtraction(x, y, z);
  printAndCompare(dyn, sta, "subtraction");
  sta = multiplication(a, b, c);
  dyn = multiplication(x, y, z);
  printAndCompare(dyn, sta, "multiplication");
  sta = division(a, b, c);
  dyn = division(x, y, z);
  printAndCompare(dyn, sta, "division");
  sta = powTest(a, b);
  dyn = powTest(x, y);
  printAndCompare(dyn, sta, "pow1");
  sta = powTest2(b);
  dyn = powTest2(y);
  printAndCompare(dyn, sta, "pow2");
  sta = powTest3(b);
  dyn = powTest3(y);
  printAndCompare(dyn, sta, "pow3");
  sta = powTest4(b);
  dyn = powTest4(y);
  printAndCompare(dyn, sta, "pow4");
  sta = complex(a, b, c);
  dyn = complex(x, y, z);
  printAndCompare(dyn, sta, "complex");
  std::cout.flags(f);
  return 0;
}

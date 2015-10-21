#include <mthAD.h>
#include <iostream>
#include <cassert>
#include <string>
using namespace mth;
template <unsigned int N>
void printAndCompare(AD<> dyn, AD<N> sta, std::string type)
{
  std::cout << "Testing output values for " << type << '\n';
  std::cout <<"Static value is:\t" << sta.val() << '\n';
  std::cout << "Dynamic value is:\t" << dyn.val() << '\n';
  assert(dyn.val() - sta.val() < .001);
  for(int i = 0; i < N; i++)
  {
    std::cout << i << " static derivate is:\t" << sta.dx(i) << '\n';
    std::cout << i << " dyniamic derivate is:\t" << dyn.dx(i) << '\n';
    assert(dyn.dx(i) - sta.dx(i) < .001);
  }
}
template <unsigned int N>
AD<N> addition(AD<N> a, AD<N> b, AD<N> c)
{
  return a + b + c + 37;
}

template <unsigned int N>
AD<N> subtraction(AD<N> a, AD<N> b, AD<N> c)
{
  return a - b - c - 37;
}

template <unsigned int N>
AD<N> multiplication(AD<N> a, AD<N> b, AD<N> c)
{
  return a * b * c * 37;
}

template <unsigned int N>
AD<N> division(AD<N> a, AD<N> b, AD<N> c)
{
  return a / b / c;
}

template <unsigned int N>
AD<N> powTest(AD<N> a, AD<N> b)
{
  return pow(a, b);
}

template <unsigned int N>
AD<N> powTest2(AD<N> a)
{
  return pow(a, 3.0);
}

template <unsigned int N>
AD<N> powTest3(AD<N> a)
{
  return pow(3.0, a);
}

template <unsigned int N>
AD<N> powTest4(AD<N> a)
{
  return pow(3, sin(a));
}

template <unsigned int N>
AD<N> complex(AD<N> a, AD<N> b, AD<N> c)
{
  return pow(a, 12) * sqrt(b / c) / sin(c * b) * 37 / c;
}

int main(int argc, char const *argv[])
{
  AD<3> a = 1.0;
  AD<3> b = 3.0;
  AD<3> c = 7.0;
  a.diff(0);
  b.diff(1);
  c.diff(2);
  AD<> x = 1.0;
  AD<> y = 3.0;
  AD<> z = 7.0;
  x.diff(0);
  y.diff(1);
  z.diff(2);
  printAndCompare(x, a, "sanity check");
  printAndCompare(y, b, "sanity check");
  printAndCompare(z, c, "sanity check");
  AD<3> sta = addition(a, b, c);
  AD<> dyn = addition(x, y, z);
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
  sta = complex(complex(a, b, c), a, c);
  dyn = complex(complex(x, y, z), x, z);
  printAndCompare(dyn, sta, "recursive complex");
  return 0;
}
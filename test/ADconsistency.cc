#include "mthAD.h"
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <string>
using namespace mth;
typedef AD<AD<AD<double, 3>, 3>, 3> AD3;
typedef AD<AD<AD<double, 0>, 0>, 0> AD3_dyn;

void Compare3rdOrder(AD3 sta, AD3_dyn dyn, std::string equation)
{
  std::cout << "\nTesting output values for " << equation << '\n';
  for(unsigned int i = 0; i < 3; i++)
  {
    assert(std::abs((double)dyn.dx(i) - (double)sta.dx(i)) < 1e-14);
    for(unsigned int j = 0; j < 3; j++)
    {
      assert(std::abs((double)dyn.dx(i).dx(j) - (double)sta.dx(i).dx(j)) < 1e-14);
      for(unsigned int k = 0; k < 3; k++)
      {
        assert(std::abs((double)dyn.dx(i).dx(j).dx(k) 
                      - (double)sta.dx(i).dx(j).dx(k)) < 1e-14);
      }
    }
  }
}
void printAndCompare(AD<double, 0> dyn, AD<double, 3> sta, std::string type)
{
  std::cout << "Testing output values for " << type << '\n';
  std::cout << "Static value is: \t\t" << sta.val() << '\n';
  std::cout << "Dynamic value is: \t\t" << dyn.val() << '\n';
  assert(std::abs(dyn.val() - sta.val()) < 1e-14);
  for(unsigned int i = 0; i < 3; i++)
  {
    std::cout << i << " static derivative is: \t" << sta.dx(i) << '\n';
    std::cout << i << " dynamic derivative is:\t" << dyn.dx(i) << '\n';
    assert(std::abs((double)dyn.dx(i) - (double)sta.dx(i)) < 1e-14);
  }
}

template <unsigned int N>
void printThirdOrder(AD<AD<AD<double, N>, N>, N> x, std::string equation)
{
  std::cout << "\n\nEquation: " << equation << '\n';
  std::cout << "Value: " << (double)x.val() << '\n';
  std::cout << "dx: " << (double)x.dx(0) << '\n';
  std::cout << "dy: " << (double)x.dx(1) << '\n';
  std::cout << "dz: " << (double)x.dx(2) << '\n';
  std::cout << "dxx: " << (double)x.dx(0).dx(0) << '\n';
  std::cout << "dxy: " << (double)x.dx(0).dx(1) << '\n';
  std::cout << "dxz: " << (double)x.dx(0).dx(2) << '\n';
  std::cout << "dyy: " <<  (double)x.dx(1).dx(1) << '\n';
  std::cout << "dyz: " << (double)x.dx(1).dx(2) << '\n';
  std::cout << "dzz: " << (double)x.dx(2).dx(2) << '\n';
  std::cout << "dxxx: " << (double)x.dx(0).dx(0).dx(0) << '\n';
  std::cout << "dyyy: " << (double)x.dx(1).dx(1).dx(1) << '\n';
  std::cout << "dzzz: " << (double)x.dx(2).dx(2).dx(2) << '\n';
  std::cout << "dxyz: " << (double)x.dx(0).dx(1).dx(2) << "\n";

}

template <unsigned int N>
void print4thVariable(AD<AD<AD<double, N>, N>, N> x, std::string equation)
{
  printThirdOrder(x, equation);
  std::cout << "dw: " << (double)x.dx(3) << '\n';
  std::cout << "dxw: " << (double)x.dx(0).dx(3) << '\n';
  std::cout << "dyw: " << (double)x.dx(1).dx(3) << '\n';
  std::cout << "dzw: " << (double)x.dx(2).dx(3) << '\n';
  std::cout << "dww: " << (double)x.dx(3).dx(3) << '\n';
  std::cout << "dwww: " << (double)x.dx(3).dx(3).dx(3) << '\n';
}

template <class T, unsigned int N>
AD<T, N> addition(AD<T, N> a, AD<T, N> b, AD<T, N> c)
{
  return a + b + c + AD<T, N>(37);
}

template <class T, unsigned int N>
AD<T, N> subtraction(AD<T, N> a, AD<T, N> b, AD<T, N> c)
{
  return a - b - c;
}

template <class T, unsigned int N>
AD<T, N> multiplication(AD<T, N> a, AD<T, N> b, AD<T, N> c)
{
  return a * b * c;
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
  return pow(a, 12) * sqrt(b / c) / sin(c * b) / c;
}


int main()
{
  /*
  Dynamic AD consistency tests
  */
  std::ios::fmtflags f( std::cout.flags() );
  std::cout << std::setprecision(20);
  AD<double, 3> a = AD<double, 3>(1.0);
  AD<double, 3> b = AD<double, 3>(3.0);
  AD<double, 3> c = AD<double, 3>(7.0);
  a.diff(0);
  b.diff(1);
  c.diff(2);
  AD<double, 0> x = AD<double, 0>(1.0);
  AD<double, 0> y = AD<double, 0>(3.0);
  AD<double, 0> z = AD<double, 0>(7.0);
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

  /*
  Third Order Tests
  */
  AD3 X = AD3(3);
  AD3 Y = AD3(7);
  AD3 Z = AD3(5);
  X.diff(0);
  Y.diff(1);
  Z.diff(2);
  printThirdOrder(X*X*X, "X*X*X");
  printThirdOrder(X - X - X, "X-X-X");
  printThirdOrder(X - Y - Z, "X-Y-Z");
  printThirdOrder(Z*(X-Y), "Z*(X - Y)");
  printThirdOrder(X*Y*Z, "XYZ");
  printThirdOrder(AD3(35)*X*Y*Y*Y + Z*X, "35*X*Y^3 + XZ");
  printThirdOrder(X/Y, "X/Y");
  printThirdOrder(X/Z/Y, "X/Z/Y");
  printThirdOrder(exp(X), "exp(X)");
  printThirdOrder(exp(X*Y/AD3(7)), "exp(XY)");
  printThirdOrder(sin(Z), "sin(Z)");
  printThirdOrder(sin(X*Y*Z/AD3(15)), "sin(XYZ/15)");
  printThirdOrder(cos(Y), "cos(Y)");
  printThirdOrder(pow(X, 3), "X^3");
  printThirdOrder(pow(X, .5), "X^.5");
  printThirdOrder(sqrt(X), "sqrt(X)");
  printThirdOrder(tan(X), "tan(X)");
  printThirdOrder(1./sin(X), "1/sin(X)");
  printThirdOrder(1./sqrt(X), "1/sqrt(X)");
  printThirdOrder(log(X), "log(X)");
  printThirdOrder(complex(X, Y, Z), "Complex Test");

  /*
  Dynamic Third Order Tests
  */
  AD3_dyn A = AD3_dyn(3);
  AD3_dyn B = AD3_dyn(7);
  AD3_dyn C = AD3_dyn(5);
  A.diff(0, 3);
  B.diff(1, 3);
  C.diff(2, 3);
  printThirdOrder(A*A*A, "A*A*A Dynamic");
  printThirdOrder(A - A - A, "A-A-A Dynamic");
  printThirdOrder(A - B - C, "A-B-C Dynamic");
  printThirdOrder(C*(A-B), "C*(A - B) Dynamic");
  printThirdOrder(A*B*C, "ABC Dynamic");
  printThirdOrder(AD3_dyn(35)*A*B*B*B + C*A, "35*A*B^3 + AC Dynamic");
  printThirdOrder(A/B, "A/B Dynamic");
  printThirdOrder(A/C/B, "A/C/B Dynamic");
  printThirdOrder(exp(A), "exp(A) Dynamic");
  printThirdOrder(exp(A*B/AD3_dyn(7)), "exp(AB/7) Dynamic");
  printThirdOrder(sin(C), "sin(C) Dynamic");
  printThirdOrder(sin(A*B*C/AD3_dyn(15)), "sin(ABC/15) Dynamic");
  printThirdOrder(cos(B), "cos(B) Dynamic");
  printThirdOrder(pow(A, 3), "A^3 Dynamic");
  printThirdOrder(pow(A, .5), "A^.5 Dynamic");
  printThirdOrder(sqrt(A), "sqrt(A) Dynamic");
  printThirdOrder(tan(A), "tan(A) Dynamic");
  printThirdOrder(1./sin(A), "1/sin(A) Dynamic");
  printThirdOrder(1./sqrt(A), "1/sqrt(A) Dynamic");  
  printThirdOrder(log(A), "log(A)");
  printThirdOrder(complex(A, B, C), "Complex Test Dynamic");

  // Add another variable to the system.
  AD3_dyn D = AD3_dyn(4);
  D.diff(3, 4);
  print4thVariable(A * D, "A*D");
  print4thVariable(B/D, "B/D");
  print4thVariable(pow(C, D), "C^D");
  print4thVariable(pow(D, 4), "D^4");
  print4thVariable(complex(A, B, D), "Complex");

  //Testing 3rd order consistency
  Compare3rdOrder(X, A, "X");
  Compare3rdOrder(X/Y, A/B, "X/Y");
  Compare3rdOrder(tan(Z), tan(C), "tan(Z)");
  Compare3rdOrder(exp(Y), exp(B), "exp(Y)");
  Compare3rdOrder(1./sqrt(X), 1./sqrt(A), "1/sqrt(X)");
  Compare3rdOrder(complex(X, Y, Z), complex(A, B, C), "Complex");
  Compare3rdOrder(log(Z*X*Y), log(C*A*B), "log(Z*X*Y)");

  //test third order derivatives to knows outputs
  std::cout <<"\n\nTesting against Known Values\n";
  std::cout << "d/dx of X*X*X" << '\n';
  assert(std::abs(double((X*X*X).dx(0)) - 27) < 1e-14);
  std::cout << "d^3/dy of Y^2 + Y^3" << '\n';
  assert(std::abs(double((Y*Y + Y*Y*Y).dx(1).dx(1).dx(1)) - 6) < 1e-14);
  std::cout << "d^2/(dx^2) of sin(X)" << '\n';
  assert(std::abs(double((sin(X)).dx(0).dx(0)) - (-std::sin(3))) < 1e-14);
  std::cout << "d^2/(dxdy) of log(XY)" << '\n';
  assert(std::abs(double((log(X*Y)).dx(0).dx(1)) - 0) < 1e-14);
  std::cout << "d^3/(dxdydz) of exp(XYZ)" << '\n';
  assert(std::abs(double(
        (exp(X*Y*Z)).dx(0).dx(1).dx(2)) - 11341*std::exp(105)) < 1e-14);
  std::cout << "d^2/dxdy of pow(XY, .3)" << '\n';
  assert(std::abs(double(pow(X*Y, .3).dx(0).dx(1)) 
                  - 0.09/std::pow(21, .7)) < 1e-14);
  std::cout <<"d^3/dz^3 of 1/Z" << '\n';
  assert(std::abs(double((1./Z).dx(2).dx(2).dx(2)) - (-0.0096)) < 1e-14);
  return 0;
}

#include "mth_def.h"
#include <iostream>

namespace {

static void setTensor2(mth::Tensor<double>& a)
{
  a(0,0) = 2.0; a(0,1) = -1.0;
  a(1,0) = 1.0; a(1,1) = 2.0;
  std::cout << "2x2 matrix" << std::endl;
  std::cout << a << std::endl;
}

static void setTensor3(mth::Tensor<double>& a)
{
  a(0,0) = 2.0; a(0,1) = -1.0; a(0,2) = 0.0;
  a(1,0) = 1.0; a(1,1) = 2.0; a(1,2) = -1.0;
  a(2,0) = 0.0; a(2,1) = -1.0; a(2,2) = 3.0;
  std::cout << "3x3 matrix" << std::endl;
  std::cout << a << std::endl;
}

static void testMultiply(mth::Tensor<double>& a)
{
  mth::Tensor<double> b(2);
  b = a*a;
  std::cout << "\nmultiply:" << std::endl;
  std::cout << b << std::endl;
}

static void testNorm(mth::Tensor<double>& a)
{
  double n = mth::norm(a);
  std::cout << "\nnorm:" << std::endl;
  std::cout << n << std::endl;
}

static void testDeterminant(mth::Tensor<double>& a)
{
  double det = mth::determinant(a);
  std::cout << "\ndeterminant: " << std::endl;
  std::cout << det << std::endl;
}

static void testTranspose(mth::Tensor<double>& a)
{
  mth::Tensor<double> b;
  mth::transpose(a, b);
  std::cout << "\ntranspose:" << std::endl;
  std::cout << b << std::endl;
}

static void testInverse(mth::Tensor<double>& a)
{
  mth::Tensor<double> b;
  mth::inverse(a, b);
  std::cout << "\ninverse:" << std::endl;
  std::cout << b << std::endl;
}

static void testDeviatoric(mth::Tensor<double>& a)
{
  mth::Tensor<double> b;
  mth::deviatoric(a, b);
  std::cout << "\ndeviatoric:" << std::endl;
  std::cout << b << std::endl;
}

}

int main()
{
  mth::Tensor<double> a(2);
  mth::Tensor<double> b(3);
  setTensor2(a);
  setTensor3(b);
  testMultiply(a);
  testMultiply(b);
  testNorm(a);
  testNorm(b);
  testDeterminant(a);
  testDeterminant(b);
  testTranspose(a);
  testTranspose(b);
  testInverse(a);
  testInverse(b);
  testDeviatoric(a);
  testDeviatoric(b);
}

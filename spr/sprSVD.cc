/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "spr.h"

namespace spr {

int getMax(int a, int b)
{
  double r = (a > b) ? a : b;
  return r;
}

double getSign(double x)
{
  double r = (x >= 0.0) ? 1.0 : -1.0;
  return r;
}

void identitize(apf::DynamicMatrix& I)
{
  for (std::size_t i=0; i < I.getRows(); ++i)
  {
    for (std::size_t j=0; j < I.getColumns(); ++j)
      I(i,j) = 0.0;
    I(i,i) = 1.0;
  }
}

/** @brief singular value decomposition of mxn matrix A
  * @param A (In) mxn matrix to be decomposed
  *          (Out) mxn truncated orthogonal matrix U
  * @param W (Out) nx1 vector of singular values 
  * @param V (Out) nxn orthogonal matrix V 
  * @details Reduced singular value decomposition of 
  *          matrix A = U*W*V^T found via one-sided
  *          Jacobi iteration. See Algorithm 4.1 of
  *          "Jacobi's method is more accurate than QR" 
  *          by Demmel and Vesselic
  */
void decomposeSVD(apf::DynamicMatrix& A,
                  apf::DynamicVector& W,
                  apf::DynamicMatrix& V) 
{
  std::size_t m = A.getRows();
  std::size_t n = A.getColumns();
  V.setSize(n,n);
  identitize(V);
  W.setSize(n);
  int count = 1;
  int sweep = 0;
  int sweep_max = 5*int(n);
  sweep_max = getMax(sweep_max,12);
  double eps = 2.2204460492503131e-16;
  double tol = 10*double(m)*eps;
  while (count > 0 && sweep <= sweep_max)
  {
    count = n*(n-1)/2;
    for (std::size_t j=0; j < n-1; ++j)
    {
      for (std::size_t k=j+1; k < n; ++k)
      {
        double a=0, b=0, c=0;
        for (std::size_t i=0; i < m; ++i)
        {
          a += A(i,j)*A(i,j);
          b += A(i,k)*A(i,k);
          c += A(i,j)*A(i,k);
        }
        if (fabs(c)/sqrt(a*b) <= tol)
        {
          count--;
          continue;
        }
        double e = (b-a)/(2*c);
        double t = getSign(e)/(fabs(e) + sqrt(1 + e*e));
        double cs = 1/sqrt(1 + t*t);
        double sn = cs*t;
        for (std::size_t i=0; i < m; ++i)
        {
          double aij = A(i,j);
          double aik = A(i,k);
          A(i,j) = cs*aij - sn*aik;
          A(i,k) = sn*aij + cs*aik;
        }
        for (std::size_t i=0; i < n; ++i)
        {
          double vij = V(i,j);
          double vik = V(i,k);
          V(i,j) = cs*vij - sn*vik;
          V(i,k) = sn*vij + cs*vik;
        }
      }
    }
    sweep++;
  }
  if (count > 0)
    apf::fail("SPR: Jacobi iterations did not reach tolerance");
  apf::DynamicVector aj;
  for (std::size_t j=0; j < n; ++j)
  {
    A.getColumn(j,aj);
    double norm = aj.getLength();
    W(j) = norm;
    if (norm <= tol)
      aj.zero();
    else
      aj /= norm;
    A.setColumn(j,aj);
  }
}

/** @brief solve a diagonal linear system Ax = b
  * @param A (In) nxn diagonal matrix stored in vector
  * @param x (Out) nx1 solution vector
  * @param b (In) nx1 right hand side
  */
void solveDiag(apf::DynamicVector& A,
               apf::DynamicVector& x,
               apf::DynamicVector& b)
{
  std::size_t n = A.getSize();
  x.setSize(n);
  double tol = 1.0e-9;
  for (std::size_t i=0; i < n; ++i)
  {
    if (A(i) < tol)
      x(i) = 0.0;
    else
      x(i) = b(i)/A(i);
  }
}

void solveSVD(apf::DynamicMatrix& A,
              apf::DynamicVector& x,
              apf::DynamicVector& b)
{
  apf::DynamicVector W;
  apf::DynamicMatrix V;
  decomposeSVD(A,W,V);
  apf::DynamicMatrix UT;
  apf::transpose(A,UT);
  apf::DynamicVector UTb;
  apf::multiply(UT,b,UTb);
  apf::DynamicVector y;
  solveDiag(W,y,UTb);
  apf::multiply(V,y,x);
}

}

/******************************************************************************

  Copyright 2015 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#ifndef MTH_DEF_H
#define MTH_DEF_H

#include "mth.h"

namespace mth {

template <class T, unsigned M>
T norm(Vector<T,M> const& a)
{
  return sqrt(a * a);
}

template <class T>
Vector<T,3> cross(Vector<T,3> const& a, Vector<T,3> const& b)
{
  Vector3<T> r;
  r(0) = a(1)*b(2) - a(2)*b(1);
  r(1) = a(2)*b(0) - a(0)*b(2);
  r(2) = a(0)*b(1) - a(1)*b(0);
  return r;
}

template <class T>
Matrix<T,3,3> cross(Vector<T,3> const& a)
{
  Matrix3x3<T> r(
      0, -a(2),  a(1),
   a(2),     0, -a(0),
  -a(1),  a(0),     0);
  return r;
}

template <class T, unsigned N>
Vector<T,N> project(Vector<T,N> const& a, Vector<T,N> const& b)
{
  return b*((a*b)/(b*b));
}

template <class T, unsigned N>
Vector<T,N> reject(Vector<T,N> const& a, Vector<T,N> const& b)
{
  return a - project(a, b);
}

template <class T, unsigned M, unsigned N>
void transpose(Matrix<T,M,N> const& a,
    Matrix<T,N,M>& b)
{
  unsigned m = a.rows();
  unsigned n = a.rows();
  b.resize(m, n);
  for (unsigned i=0; i < m; ++i)
  for (unsigned j=0; j < n; ++j)
    b(j,i) = a(i,j);
}

template <class T>
T determinant(Matrix<T,2,2> const& a)
{
  return a(0,0)*a(1,1) - a(1,0)*a(0,1);
}

template <class T>
T determinant(Matrix<T,3,3> const& a)
{
  return
    a(0,0) * (a(1,1)*a(2,2) - a(2,1)*a(1,2)) -
    a(0,1) * (a(1,0)*a(2,2) - a(2,0)*a(1,2)) +
    a(0,2) * (a(1,0)*a(2,1) - a(2,0)*a(1,1));
}

template <class T>
Matrix<T,2,2> inverse(Matrix<T,2,2> const& a)
{
  Matrix<T,2,2> r;
  r(0,0) =  a(1,1); r(0,1) = -a(0,1);
  r(1,0) = -a(1,0); r(1,1) =  a(0,0);
  return r / determinant(a);
}

template <class T>
Matrix<T,3,3> inverse(Matrix<T,3,3> const& a)
{
  Matrix<T,3,3> r;
  Matrix<T,3,3> x = transpose(a);
  r[0] = cross(x[1], x[2]);
  r[1] = cross(x[2], x[0]);
  r[2] = cross(x[0], x[1]);
  return r / determinant(a);
}

template <class T>
T trace(Tensor<T> const& a)
{
  unsigned dim = a.dim();
  T t = a(0,0);
  for (unsigned i=1; i < dim; ++i)
    t += a(i,i);
  return t;
}

template <class T, unsigned M, unsigned N>
void multiply(Matrix<T,M,N> const& a, Vector<T,N> const& b,
    Vector<T,M>& c)
{
  unsigned m = a.rows();
  unsigned n = a.cols();
  c.resize(m);
  for (unsigned i = 0; i < m; ++i) {
    c(i) = 0;
    for (unsigned j = 0; j < n; ++j)
      c(i) += a(i,j) * b(j);
  }
}

template <class T, unsigned M, unsigned N, unsigned O>
void multiply(Matrix<T,M,O> const& a, Matrix<T,O,N> const& b,
    Matrix<T,M,N>& c)
{
  unsigned m = a.rows();
  unsigned n = b.cols();
  unsigned o = a.cols();
  c.resize(m, n);
  for (unsigned i = 0; i < m; ++i)
  for (unsigned j = 0; j < n; ++j) {
    c(i,j) = 0;
    for (unsigned l = 0; l < o; ++l)
      c(i,j) += a(i,l) * b(l,j);
  }
}

}

#endif

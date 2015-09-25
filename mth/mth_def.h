/******************************************************************************

  Copyright 2015 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  This work is open source software, licensed under the terms of the
  BSD license as described in the LICENSE file in the top-level directory.

*******************************************************************************/

#ifndef MTH_DEF_H
#define MTH_DEF_H

#include <cassert>

namespace mth {

template <class T>
Vector3<T> cross(Vector3<T> const& a, Vector3<T> const& b)
{
  Vector3<T> r;
  r[0] = a[1]*b[2] - a[2]*b[1];
  r[1] = a[2]*b[0] - a[0]*b[2];
  r[2] = a[0]*b[1] - a[1]*b[0];
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

template <class T, unsigned M>
Matrix<T,M,M> eye(unsigned m)
{
  if (M != 0) assert(m == M);
  Matrix<T,M,M> r(m,m);
  r.zero();
  for (unsigned i=0; i < m; ++i)
    r(i,i) = (T)1.0;
  return r;
}

template <class T, unsigned M, unsigned N>
Matrix<T,M,N> transpose(Matrix<T,M,N> const& a)
{
  unsigned m = a.rows();
  unsigned n = a.cols();
  Matrix<T,M,N> r(m,n);
  for (unsigned i=0; i < m; ++i)
  for (unsigned j=0; j < n; ++j)
    r(j,i) = a(i,j);
  return r;
}

template <class T, unsigned M>
T trace(Matrix<T,M,M> const& a)
{
  unsigned m = a.rows();
  unsigned n = a.cols();
  assert(m == n);
  T t = a(0,0);
  for (unsigned i=1; i < m; ++i)
    t += a(i,i);
  return t;
}

template <class T, unsigned M>
T determinant(Matrix<T,M,M> const& a)
{
  unsigned m = a.rows();
  unsigned n = a.cols();
  assert(m == n);
  T d = (T)0.0;
  switch (m)
  {
    case 2:
       d = a(0,0)*a(1,1) - a(0,1)*a(1,0);
       break;
    case 3:
       d = 
         a(0,0)*(a(1,1)*a(2,2) - a(1,2)*a(2,1)) -
         a(0,1)*(a(1,0)*a(2,2) - a(1,2)*a(2,0)) +
         a(0,2)*(a(1,0)*a(2,1) - a(1,1)*a(2,0));
       break;
    default:
      fprintf(stderr,"det: unsupported dim\n");
      abort();
  }
  return d;
}

template <class T, unsigned M>
Matrix<T,M,M> inverse(Matrix<T,M,M> const& a)
{
  unsigned m = a.rows();
  unsigned n = a.cols();
  assert(m == n);
  Matrix<T,M,M> r(m,m);
  switch (m)
  {
    case 2:
      r(0,0) =  a(1,1);  r(0,1) = -a(0,1);
      r(1,0) = -a(1,0);  r(1,1) =  a(0,0);
      r /= det(a);
      break;
    default:
      fprintf(stderr,"inverse: unsupported dim\n");
      abort();
  }
  return r;
}

template <class T, unsigned M>
Matrix<T,M,M> deviatoric(Matrix<T,M,M> const& a)
{
  unsigned m = a.rows();
  unsigned n = a.cols();
  assert(m == n);
  Matrix<T,M,M> r(m,m);
  r = a;
  T t = trace(a) / (T)m;
  for (unsigned i=0; i < m; ++i)
    r(i,i) -= t;
  return r;
}

template <class T, unsigned M, unsigned N>
T norm(Matrix<T,M,N> const& a)
{
  unsigned m = a.rows();
  unsigned n = a.cols();
  T r = (T)0.0;
  for (unsigned i=0; i < m; ++i)
  for (unsigned j=0; j < n; ++j)
    r += a(i,j)*a(i,j);
  return sqrt(r);
}

}

#endif

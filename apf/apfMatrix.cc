/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfMatrix.h"
#include "apf2mth.h"
#include <mthQR.h>
#include <cassert>

namespace apf {

Matrix3x3 cross(Vector3 const& u)
{
  return Matrix3x3( 0    ,-u.z(), u.y(),
                    u.z(), 0    ,-u.x(),
                   -u.y(), u.x(), 0    );
}

Matrix3x3 rotate(Vector3 const& u, double a)
{
  Matrix3x3 I(1,0,0,
              0,1,0,
              0,0,1);
  return I*cos(a) + cross(u)*sin(a) + tensorProduct(u,u)*(1-cos(a));
}

/* this is defined outside getFrame because in
   C++98 local types cannot be template arguments,
   and we use the templated std::swap function */
struct SortStruct
{
  int i;
  double m;
};

Matrix3x3 getFrame(Vector3 const& v)
{
  Matrix<3,3> A;
  A[0] = v;
  /* tiny custom code to sort components by absolute value */
  SortStruct s[3] =
  {{0,fabs(v[0])},{1,fabs(v[1])},{2,fabs(v[2])}};
  if (s[2].m > s[1].m)
    std::swap(s[1],s[2]);
  if (s[1].m > s[0].m)
    std::swap(s[0],s[1]);
  if (s[2].m > s[1].m)
    std::swap(s[1],s[2]);
  /* done, components sorted by increasing magnitude */
  int lc = s[0].i;
  int mc = s[1].i;
  int sc = s[2].i;
  /* use the 2D rotation on the largest components
     (rotate v around the smallest axis) */
  A[1][lc] = -v[mc];
  A[1][mc] = v[lc];
  /* and make the last component zero so that A[0] * A[1] == 0 */
  A[1][sc] = 0;
  /* now we have 2 orthogonal (though not unit) vectors, cross
     product gives the third */
  A[2] = cross(A[0],A[1]);
  return A;
}

int eigen(Matrix3x3 const& A,
          Vector<3>* eigenVectors,
          double* eigenValues)
{
  mth::Matrix<double,3,3> A2 = to_mth(A);
  mth::Matrix<double,3,3> L;
  mth::Matrix<double,3,3> Q;
  bool converged = mth::eigenQR(A2, L, Q, 100);
  assert(converged);
  for (unsigned i = 0; i < 3; ++i)
    eigenValues[i] = L(i,i);
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j)
    eigenVectors[j][i] = Q(i,j);
  return 3;
}

template <std::size_t M, std::size_t N>
Matrix<M - 1, N - 1> getMinor(Matrix<M,N> const& A,
    std::size_t i, std::size_t j)
{
  Matrix<N - 1, M - 1> B;
  std::size_t m = 0;
  for (std::size_t k = 0; k < M; ++k)
    if (k != i) {
      std::size_t n = 0;
      for (std::size_t l = 0; l < N; ++l)
        if (l != j) {
          B[m][n] = A[k][l];
          ++n;
        }
      ++m;
    }
  return B;
}

template <std::size_t M, std::size_t N>
double getCofactor(Matrix<M,N> const& A, std::size_t i, std::size_t j)
{
  Matrix<M - 1, N - 1> B = getMinor(A, i, j);
  double dM = getDeterminant(B);
  double sign = ((i + j) % 2) ? -1 : 1;
  return sign * dM;
}

template <std::size_t M, std::size_t N>
double getDeterminant(Matrix<M,N> const& A)
{
  double d = 0;
  for (std::size_t i = 0; i < M; ++i)
    d += A[i][0] * getCofactor(A, i, 0);
  return d;
}

double getDeterminant(Matrix<1,1> const& A)
{
  return A[0][0];
}

template Matrix<1,1> getMinor(Matrix<2,2> const& A, std::size_t i, std::size_t j);
template Matrix<2,2> getMinor(Matrix<3,3> const& A, std::size_t i, std::size_t j);
template Matrix<3,3> getMinor(Matrix<4,4> const& A, std::size_t i, std::size_t j);

template double getDeterminant(Matrix<2,2> const& A);
template double getDeterminant(Matrix<3,3> const& A);
template double getDeterminant(Matrix<4,4> const& A);

}

/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfMatrix.h"
#include <complex>

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

struct Cubic
{
/* coefficients of a cubic polynomial */
  double a;
  double b;
  double c;
  double d;

  Cubic(double A, double B, double C, double D)
  {
    a = A; b = B; c = C; d = D;
  }
  Cubic() {}

/* http://en.wikipedia.org/wiki/Cubic_function#The_nature_of_the_roots */
  double getDiscriminant()
  {
    return 18 * a * b * c * d
         -  4 * b * b * b * d
         +      b * b * c * c
         -  4 * a * c * c * c
         - 27 * a * a * d * d;
  }

  struct Roots
  {
    int n;
    double x[3];
  };

/* http://en.wikipedia.org/wiki/Cubic_function#General_formula_for_roots */
/*   zeroes cause havoc, refer to the following: */
/* http://en.wikipedia.org/wiki/Cubic_function#Special_cases */
  void solve(Roots& r)
  {
    double disc = getDiscriminant();
    if (disc > 0)
      r.n = 3;
    else
      r.n = 1;
    std::complex<double> u[3];
    u[0] = 1;
    u[1] = std::complex<double>(-1,  std::sqrt(3.)) / 2.;
    u[2] = std::complex<double>(-1, -std::sqrt(3.)) / 2.;
    double d0 = b * b - 3 * a * c;
    double d1 = 2 * b * b * b
             -  9 * a * b * c
             + 27 * a * a * d;
    std::complex<double> tmp;
    if (d0) {
      tmp = d1 * d1 - 4 * d0 * d0 * d0;
      tmp = sqrt(tmp);
    } else {
      tmp = d1;
    }
    std::complex<double> C = pow((d1 + tmp) / 2., 1. / 3.);
    std::complex<double> x[3];
    if (abs(C)) {
      for (int k = 0; k < 3; ++k)
        x[k] = -(1 / (3 * a)) *
          (b + u[k] * C + (d0 / (u[k] * C)));
    } else {
      if (d0) {
        x[0] = x[1] = (9 * a * d - b * c) / (2 * d0);
        x[2] = (4 * a * b * c - 9 * a * a * d - b * b * b) / (a * d0);
      } else {
        for (int k = 0; k < 3; ++k)
          x[k] = -(b / (3 * a));
      }
    }
    if (r.n == 3)
      for (int k = 0; k < 3; ++k)
        r.x[k] = x[k].real();
    else {
      int best = 0;
      for (int k = 1; k < 3; ++k)
        if (fabs(x[k].imag()) < fabs(x[best].imag()))
          best = k;
      r.x[0] = x[best].real();
    }
  }
};

static double tr(Matrix3x3 const& m)
{
  return m[0][0] + m[1][1] + m[2][2];
}

/* http://mathworld.wolfram.com/CharacteristicPolynomial.html */
static void getCharacteristicPolynomial(Matrix3x3 const& A, Cubic& p)
{
  double tA = tr(A);
  double c2 = (1. / 2.) * (tA * tA - tr(A * A));
  p.a = -1;
  p.b = tA;
  p.c = -c2;
  p.d = getDeterminant(A);
}

static int getEigenvalues(Matrix3x3 const& A, double* v)
{
  Cubic p;
  getCharacteristicPolynomial(A, p);
  Cubic::Roots l;
  p.solve(l);
  for (int i = 0; i < l.n; ++i)
    v[i] = l.x[i];
  return l.n;
}

static void getEigenvector(Matrix3x3 const& A, double l, Vector<3>& v)
{
  Matrix3x3 eye(1,0,0,
                0,1,0,
                0,0,1);
  Matrix3x3 basis = transpose(A - eye * l);
/* the rows of this ^ matrix should in the general case
   span a plane, or in bad cases a line or nothing.
   We assume the good case and take the most outstanding cross product
   of any pair of rows */
  Vector3 c[3];
  c[0] = cross(basis[0],basis[1]);
  c[1] = cross(basis[1],basis[2]);
  c[2] = cross(basis[2],basis[0]);
  int best = 0;
  for (int i = 1; i < 3; ++i)
    if (c[i].getLength() > c[best].getLength())
      best = i;
/* if the null space is higher than 1-dimensional, this
   will break in div by zero, or divide by something very close to 0 */
  v = c[best].normalize();
}

int eigen(Matrix3x3 const& A,
          Vector<3>* eigenVectors,
          double* eigenValues)
{
  int n = getEigenvalues(A, eigenValues);
  for (int i = 0; i < n; ++i)
    getEigenvector(A, eigenValues[i], eigenVectors[i]);
  return n;
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

template <std::size_t M, std::size_t N>
std::ostream& print(std::ostream& s, Matrix<N,M> const& A)
{
  for (std::size_t i = 0; i < M; ++i) {
    for (std::size_t j = 0; j < N; ++j)
      s << A[i][j] << ' ';
    s << '\n';
  }
  return s;
}

}

std::ostream& operator<<(std::ostream& s, apf::Matrix<1,1> const& A)
{
  return apf::print(s, A);
}

std::ostream& operator<<(std::ostream& s, apf::Matrix<2,2> const& A)
{
  return apf::print(s, A);
}

std::ostream& operator<<(std::ostream& s, apf::Matrix<3,3> const& A)
{
  return apf::print(s, A);
}

std::ostream& operator<<(std::ostream& s, apf::Matrix<4,4> const& A)
{
  return apf::print(s, A);
}

#include "mthQR.h"
#include "mth_def.h"
#include <cassert>

namespace mth {

static double sign(double x)
{
  return (x < 0) ? -1 : 1;
}

template <class T>
static T square(T x)
{
  return x * x;
}

template <class T, unsigned M, unsigned N>
static bool get_reflector(
    Matrix<T,M,N> const& a,
    Vector<T,M>& v,
    unsigned k)
{
  unsigned m = a.rows();
  assert(v.size() == m);
  double cnorm = 0;
  for (unsigned i = k; i < m; ++i)
    cnorm += square(a(i,k));
  cnorm = sqrt(cnorm);
  if (cnorm < 1e-10)
    return false;
  for (unsigned i = 0; i < k; ++i)
    v(i) = 0;
  for (unsigned i = k; i < m; ++i)
    v(i) = a(i,k);
  v(k) += sign(a(k,k)) * cnorm;
  double rnorm = 0;
  for (unsigned i = k; i < m; ++i)
    rnorm += square(v(i));
  rnorm = sqrt(rnorm);
  for (unsigned i = k; i < m; ++i)
    v(i) /= rnorm;
  return true;
}

template <class T, unsigned M, unsigned N>
static void reflect_columns(
    Vector<T,M> const& v,
    Matrix<T,M,N>& a,
    unsigned k)
{
  unsigned m = a.rows();
  unsigned n = a.cols();
  for (unsigned j = k; j < n; ++j) {
    double dot = 0;
    for (unsigned i = k; i < m; ++i)
      dot += a(i,j) * v(i);
    for (unsigned i = k; i < m; ++i)
      a(i,j) -= 2 * dot * v(i);
  }
}

template <class T, unsigned M>
static void reflect_rows(
    Vector<T,M> const& v,
    Matrix<T,M,M>& q,
    unsigned k)
{
  unsigned m = q.rows();
  for (unsigned i = 0; i < m; ++i) {
    double dot = 0;
    for (unsigned j = k; j < m; ++j)
      dot += q(i,j) * v(j);
    for (unsigned j = k; j < m; ++j)
      q(i,j) -= 2 * dot * v(j);
  }
}

template <class T, unsigned M, unsigned N>
static bool qr_step(
    Matrix<T,M,N>& a,
    Matrix<T,M,M>& q,
    Vector<T,M>& v,
    unsigned k)
{
  if (!get_reflector(a, v, k))
    return false;
  reflect_columns(v, a, k);
  reflect_rows(v, q, k);
  return true;
}

template <class T, unsigned M>
static void fill_identity(Matrix<T,M,M>& q)
{
  unsigned m = q.rows();
  for (unsigned i = 0; i < m; ++i)
  for (unsigned j = 0; j < m; ++j)
    q(i,j) = ((double)(i==j));
}

template <class T, unsigned M, unsigned N>
unsigned decomposeQR(
    Matrix<T,M,N> const& a,
    Matrix<T,M,M>& q,
    Matrix<T,M,N>& r)
{
  unsigned m = a.rows();
  unsigned n = a.cols();
  assert(m >= n);
  q.resize(m, m);
  fill_identity(q);
  Vector<T,M> v_scratch;
  v_scratch.resize(m);
  r = a;
  unsigned rank = 0;
  for (unsigned k = 0; k < n; ++k)
    if (qr_step(r, q, v_scratch, k))
      ++rank;
  return rank;
}

template unsigned decomposeQR(Matrix<double,3,3> const& a,
    Matrix<double,3,3>& q, Matrix<double,3,3>& r);
template unsigned decomposeQR(Matrix<double,0,0> const& a,
    Matrix<double,0,0>& q, Matrix<double,0,0>& r);

template <class T, unsigned M, unsigned N>
void backsubUT(
    Matrix<T,M,N> const& a,
    Vector<T,M> const& b,
    Vector<T,N>& x)
{
  unsigned m = a.rows();
  unsigned n = a.cols();
  assert(m >= n);
  x.resize(n);
  for (unsigned ii = 0; ii < n; ++ii) {
    unsigned i = n - ii - 1;
    x(i) = b(i);
    for (unsigned j = i + 1; j < n; ++j)
      x(i) -= a(i,j) * x(j);
    x(i) /= a(i,i);
  }
}

template void backsubUT(Matrix<double,0,0> const& a,
    Vector<double,0> const& b, Vector<double,0>& x);

template <class T, unsigned M, unsigned N>
void solveFromQR(Matrix<T,M,M> const& q,
    Matrix<T,M,N> const& r,
    Vector<T,M> const& b, Vector<T,N>& x)
{
  Matrix<T,M,M> qt;
  transpose(q, qt);
  Vector<T,M> y;
  multiply(qt, b, y);
  backsubUT(r, y, x);
}

template void solveFromQR(Matrix<double,0,0> const& q,
    Matrix<double,0,0> const& r,
    Vector<double,0> const& b, Vector<double,0>& x);

template <class T, unsigned M, unsigned N>
bool solveQR(Matrix<T,M,N> const& a,
    Vector<T,M> const& b, Vector<T,N>& x)
{
  Matrix<T,M,M> q;
  Matrix<T,M,N> r;
  unsigned rank = decomposeQR(a, q, r);
  if (rank != a.cols())
    return false;
  solveFromQR(q, r, b, x);
  return true;
}

template bool solveQR(Matrix<double,0,0> const& a, Vector<double,0> const& b,
    Vector<double,0>& x);

}

#include "mthQR.h"
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
  r = a;
  fill_identity(q);
  unsigned m = a.rows();
  q.resize(m, m);
  Vector<T,M> v_scratch;
  v_scratch.resize(m);
  unsigned rank = 0;
  for (unsigned k = 0; k < m; ++k)
    if (qr_step(r, q, v_scratch, k))
      ++rank;
  return rank;
}

template unsigned decomposeQR(Matrix<double,3,3> const& a,
    Matrix<double,3,3>& q, Matrix<double,3,3>& r);
template unsigned decomposeQR(Matrix<double,0,0> const& a,
    Matrix<double,0,0>& q, Matrix<double,0,0>& r);

}

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
    unsigned k,
    unsigned o = 0)
{
  unsigned m = a.rows();
  v.resize(m);
  double cnorm = 0;
  for (unsigned i = k + o; i < m; ++i)
    cnorm += square(a(i, k));
  cnorm = sqrt(cnorm);
  if (cnorm < 1e-10)
    return false;
  for (unsigned i = 0; i < k + o; ++i)
    v(i) = 0;
  for (unsigned i = k + o; i < m; ++i)
    v(i) = a(i, k);
  v(k + o) += sign(a(k + o, k)) * cnorm;
  double rnorm = 0;
  for (unsigned i = k + o; i < m; ++i)
    rnorm += square(v(i));
  rnorm = sqrt(rnorm);
  for (unsigned i = k + o; i < m; ++i)
    v(i) /= rnorm;
  return true;
}

template <class T, unsigned M, unsigned N>
static void reflect_columns(
    Vector<T,M> const& v,
    Matrix<T,M,N>& a,
    unsigned k,
    unsigned o = 0)
{
  unsigned m = a.rows();
  unsigned n = a.cols();
  for (unsigned j = 0; j < n; ++j) {
    double dot = 0;
    for (unsigned i = k + o; i < m; ++i)
      dot += a(i,j) * v(i);
    for (unsigned i = k + o; i < m; ++i)
      a(i,j) -= 2 * dot * v(i);
  }
}

template <class T, unsigned M>
static void reflect_rows(
    Vector<T,M> const& v,
    Matrix<T,M,M>& q,
    unsigned k,
    unsigned o = 0)
{
  unsigned m = q.rows();
  for (unsigned i = 0; i < m; ++i) {
    double dot = 0;
    for (unsigned j = k + o; j < m; ++j)
      dot += q(i,j) * v(j);
    for (unsigned j = k + o; j < m; ++j)
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

template <class T, unsigned M>
void reduceToHessenberg(Matrix<T,M,M> const& a, Matrix<T,M,M>& q,
    Matrix<T,M,M>& h)
{
  unsigned m = a.rows();
  assert(a.cols() == m);
  q.resize(m,m);
  fill_identity(q);
  h = a;
  Vector<T,M> v;
  for (unsigned k = 0; k < m - 2; ++k)
    if (get_reflector(h, v, k, 1)) {
      reflect_columns(v, h, k, 1);
      reflect_rows(v, h, k, 1);
      reflect_rows(v, q, k, 1);
    }
}

template void reduceToHessenberg(Matrix<double,3,3> const& a,
    Matrix<double,3,3>& q, Matrix<double,3,3>& h);

template <class T, unsigned M>
static bool reduce(Matrix<T,M,M> const& a, unsigned& red_m)
{
  while (red_m > 1) {
    if (fabs(a(red_m - 2, red_m - 1)) < 1e-10 &&
        fabs(a(red_m - 1, red_m - 2)) < 1e-10)
      --red_m;
    else
      return true;
  }
  return false;
}

template <class T, unsigned M>
static double get_wilkinson_shift(Matrix<T,M,M> const& a, unsigned red_m)
{
  double amm1 = a(red_m-2, red_m-2);
  double am = a(red_m-1, red_m-1);
  double bmm1 = a(red_m-2, red_m-1);
  double sig = (amm1 - am) / 2;
  double denom = ( fabs(sig) + sqrt(square(sig) + square(bmm1)));
  assert(fabs(denom) > 1e-10);
  double mu = am - ((sign(sig) * square(bmm1)) / denom);
  return mu;
}

template <class T, unsigned M>
static void shift_matrix(Matrix<T,M,M>& a, double mu)
{
  unsigned m = a.rows();
  for (unsigned i = 0; i < m; ++i)
    a(i,i) -= mu;
}

template <class T, unsigned M>
bool eigenQR(Matrix<T,M,M> const& a,
    Matrix<T,M,M>& l,
    Matrix<T,M,M>& q,
    unsigned max_iters)
{
  Matrix<T,M,M>& a_k = l;
  reduceToHessenberg(a, q, a_k);
  Matrix<T,M,M> r_k;
  Matrix<T,M,M> q_k;
  Matrix<T,M,M> tmp_q;
  unsigned red_m = a.rows();
  for (unsigned i = 0; i < max_iters; ++i) {
    if (!reduce(a_k, red_m))
      return true;
    double mu = get_wilkinson_shift(a_k, red_m);
    shift_matrix(a_k, mu);
    decomposeQR(a_k, q_k, r_k);
    multiply(r_k, q_k, a_k);
    shift_matrix(a_k, -mu);
    multiply(q, q_k, tmp_q);
    q = tmp_q;
  }
  return false;
}

template bool eigenQR(Matrix<double,3,3> const& a, Matrix<double,3,3>& l,
    Matrix<double,3,3>& q, unsigned max_iters);

}

#include "parma_rib.h"
#include <apfNew.h>
#include <algorithm>
#include <mthQR.h>
#include <mth_def.h>
#include <cassert>

#include <iostream>
#include <iomanip>

namespace parma {

struct Compare
{
  mth::Vector3<double> normal;
  bool operator()(Body* a, Body* b)
  {
    return (a->point * normal) < (b->point * normal);
  }
};

Bodies::Bodies(Body* arr, int n_)
{
  n = n_;
  body = new Body*[n];
  for (int i = 0; i < n; ++i)
    body[i] = arr + i;
}

Bodies::~Bodies()
{
  n = 0;
  body = NULL;
}

void Bodies::destroy()
{
  delete [] body;
}

/* http://en.wikipedia.org/wiki/Moment_of_inertia#Angular_momentum */
static mth::Matrix3x3<double> getInertiaContribution(Body const* b)
{
  mth::Matrix3x3<double> c = mth::cross(b->point);
  return c * c * -(b->mass);
}

static mth::Matrix3x3<double> getInertiaMatrix(Bodies const* b)
{
  mth::Matrix3x3<double> m(
      0,0,0,
      0,0,0,
      0,0,0);
  for (int i = 0; i < b->n; ++i)
    m = m + getInertiaContribution(b->body[i]);
  return m;
}

static mth::Matrix3x3<double> normalize(mth::Matrix3x3<double> const& A)
{
  double max = 0;
  for (unsigned i = 0; i < 3; ++i)
  for (unsigned j = 0; j < 3; ++j) {
    double val = fabs(A(i,j));
    if (val > max)
      max = val;
  }
  if (max < 1e-10)
    return A;
  return A / max;
}

static void getWeakestEigenvector(mth::Matrix3x3<double> const& A_in,
    mth::Vector3<double>& v)
{
  std::cout << std::scientific << std::setprecision(10);
  std::cerr << std::scientific << std::setprecision(10);
  mth::Matrix3x3<double> A, l, q;
  /* this will help with general conditioning of the eigenvalue solve,
     mostly superstition at this point */
  A = normalize(A_in);
  /* find the eigenvalues (l) and eigenvectors (q) of (A) */
  bool converged = mth::eigenQR(A, l, q, 100);
  assert(converged);
  unsigned best = 0;
  /* find the *smallest* eigenvalue */
  for (unsigned i = 1; i < 3; ++i)
    if (fabs(l(i,i)) < fabs(l(best,best)))
      best = i;
  /* return the corresponding eigenvector */
  for (unsigned i = 0; i < 3; ++i)
    v(i) = q(i,best);
}

static mth::Vector3<double> getBisectionNormal(Bodies const* b)
{
  mth::Matrix3x3<double> im = getInertiaMatrix(b);
  mth::Vector3<double> v;
  getWeakestEigenvector(im, v);
  return v;
}

static double getTotalMass(Bodies const* b)
{
  double mass = 0;
  for (int i = 0; i < b->n; ++i)
    mass += b->body[i]->mass;
  return mass;
}

static mth::Vector3<double> getCenterOfGravity(Bodies const* b)
{
  mth::Vector3<double> c(0,0,0);
  for (int i = 0; i < b->n; ++i) {
    Body* body = b->body[i];
    c = c + (body->point * body->mass);
  }
  return c / getTotalMass(b);
}

static void centerBodies(Bodies* b, mth::Vector3<double> const& c)
{
  for (int i = 0; i < b->n; ++i) {
    Body* body = b->body[i];
    body->point = body->point - c;
  }
}

int findSortedMedian(Bodies const* b)
{
  double total = getTotalMass(b);
  double half = 0;
  for (int i = 0; i < b->n; ++i) {
    if (half >= total / 2)
      return i;
    half += b->body[i]->mass;
  }
  return b->n;
}

void bisect(Bodies* all, Bodies* left, Bodies* right)
{
  mth::Vector3<double> c = getCenterOfGravity(all);
  centerBodies(all, c);
  Compare comp;
  comp.normal = getBisectionNormal(all);
  std::sort(all->body, all->body + all->n, comp);
  int mid = findSortedMedian(all);
  left->n = mid;
  right->n = all->n - mid;
/* in-place bisection, left and right point to the same array as all */
  left->body = all->body;
  right->body = all->body + mid;
}

void recursivelyBisect(Bodies* all, int depth, Bodies out[])
{
  if (!depth) {
    *out = *all;
    return;
  }
  Bodies left;
  Bodies right;
  bisect(all, &left, &right);
  --depth;
  recursivelyBisect(&left, depth, out);
  recursivelyBisect(&right, depth, out + (1 << depth));
}

}

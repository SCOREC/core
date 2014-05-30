#include "parma_rib.h"
#include "parma_rib_math.h"
#include <apfNew.h>
#include <algorithm>

namespace parma {

struct Compare
{
  apf::Vector3 normal;
  double median;
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

void Bodies::destroy()
{
  delete [] body;
}

/* http://en.wikipedia.org/wiki/Moment_of_inertia#Angular_momentum */
static apf::Matrix3x3 getInertiaContribution(Body const* b)
{
  apf::Matrix3x3 c = cross(b->point);
  return c * c * b->mass;
}

static apf::Matrix3x3 getInertiaMatrix(Bodies const* b)
{
  apf::Matrix3x3 m(0,0,0,
                   0,0,0,
                   0,0,0);
  for (int i = 0; i < b->n; ++i)
    m = m + getInertiaContribution(b->body[i]);
  return m * -1;
}

static apf::Vector3 getBisectionNormal(Bodies const* b)
{
  apf::Matrix3x3 im = getInertiaMatrix(b);
  apf::Vector3 v;
  getPrincipalEigenvector(im, v);
  return v;
}

static double getTotalMass(Bodies const* b)
{
  double mass = 0;
  for (int i = 0; i < b->n; ++i)
    mass += b->body[i]->mass;
  return mass;
}

static apf::Vector3 getCenterOfGravity(Bodies const* b)
{
  apf::Vector3 c(0,0,0);
  for (int i = 0; i < b->n; ++i) {
    Body* body = b->body[i];
    c = c + (body->point * body->mass);
  }
  return c / getTotalMass(b);
}

static void centerBodies(Bodies* b, apf::Vector3 const& c)
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
    half += b->body[i]->mass;
    if (half >= total / 2)
      return i;
  }
  abort();
  return -1;
}

void bisect(Bodies* all, Bodies* left, Bodies* right)
{
  apf::Vector3 c = getCenterOfGravity(all);
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

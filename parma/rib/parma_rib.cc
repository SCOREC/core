#include "parma_rib.h"
#include "parma_rib_math.h"
#include <algorithm>

namespace parma {

struct Bisector
{
  apf::Vector3 normal;
  bool operator()(Body* body)
  {
    return body->point * normal < 0;
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
  return m;
}

static apf::Vector3 getBisectionNormal(Bodies const* b)
{
  apf::Matrix3x3 im = getInertiaMatrix(b);
  apf::Vector3 v;
  getPrincipalEigenvector(im, v);
  return v;
}

static apf::Vector3 getCenterOfGravity(Bodies const* b)
{
  apf::Vector3 c(0,0,0);
  for (int i = 0; i < b->n; ++i) {
    Body* body = b->body[i];
    c = c + (body->point * body->mass);
  }
  return c / b->n;
}

static void centerBodies(Bodies* b, apf::Vector3 const& c)
{
  for (int i = 0; i < b->n; ++i) {
    Body* body = b->body[i];
    body->point = body->point - c;
  }
}

void testBisectionPlane(Body* b, int n)
{
  Bodies bodies(b, n);
  apf::Vector3 c = getCenterOfGravity(&bodies);
  std::cerr << c << '\n';
  centerBodies(&bodies, c);
  apf::Vector3 normal = getBisectionNormal(&bodies);
  std::cerr << normal << '\n';
  bodies.destroy();
}

void bisect(Bodies* all, Bodies* left, Bodies* right)
{
  apf::Vector3 c = getCenterOfGravity(all);
  centerBodies(all, c);
  Bisector b;
  b.normal = getBisectionNormal(all);
  Body** mid = std::partition(all->body, all->body + all->n, b);
  left->n = mid - all->body;
  right->n = all->n - left->n;
/* in-place bisection, left and right point to the same array as all */
  left->body = all->body;
  right->body = mid;
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

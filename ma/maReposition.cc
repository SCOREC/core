#include "maReposition.h"
#include <apfMesh.h>
#include <cmath>

#include <apf.h>
#include <sstream>

namespace ma {

/* automatic differentiation infrastructure
 * for getting the derivative of quality
 * with respect to vertex position
 */
class AD
{
  public:
    AD() {}
    AD(double v):
      v_(v)
    {
      dx_[0] = 0;
      dx_[1] = 0;
      dx_[2] = 0;
    }
    AD(double v, int dim):
      v_(v)
    {
      dx_[0] = 0;
      dx_[1] = 0;
      dx_[2] = 0;
      dx_[dim] = 1;
    }
    AD(double v, double a, double b, double c)
    {
      v_ = v;
      dx_[0] = a;
      dx_[1] = b;
      dx_[2] = c;
    }
    double val() const
    {
      return v_;
    }
    double dx(int dim) const
    {
      return dx_[dim];
    }
  private:
    double v_;
    double dx_[3];
};

static AD operator-(AD a, AD b)
{
  return AD(
      a.val() - b.val(),
      a.dx(0) - b.dx(0),
      a.dx(1) - b.dx(1),
      a.dx(2) - b.dx(2));
}

static AD operator+(AD a, AD b)
{
  return AD(
      a.val() + b.val(),
      a.dx(0) + b.dx(0),
      a.dx(1) + b.dx(1),
      a.dx(2) + b.dx(2));
}

static AD operator*(AD a, AD b)
{
  return AD(
      a.val() * b.val(),
      a.dx(0) * b.val() + a.val() * b.dx(0),
      a.dx(1) * b.val() + a.val() * b.dx(1),
      a.dx(2) * b.val() + a.val() * b.dx(2));
}

static AD operator/(AD a, AD b)
{
  return AD(
      a.val() / b.val(),
      (a.dx(0) * b.val() - a.val() * b.dx(0)) / (b.val() * b.val()),
      (a.dx(1) * b.val() - a.val() * b.dx(1)) / (b.val() * b.val()),
      (a.dx(2) * b.val() - a.val() * b.dx(2)) / (b.val() * b.val()));
}

static AD sqrt(AD a)
{
  return AD(
      std::sqrt(a.val()),
      a.dx(0) / (2. * std::sqrt(a.val())),
      a.dx(1) / (2. * std::sqrt(a.val())),
      a.dx(2) / (2. * std::sqrt(a.val())));
}

class ADVec
{
  public:
    ADVec() {}
    ADVec(AD a, AD b, AD c)
    {
      c_[0] = a;
      c_[1] = b;
      c_[2] = c;
    }
    AD operator[](int i) const
    {
      return c_[i];
    }
  private:
    AD c_[3];
};

static ADVec operator-(ADVec a, ADVec b)
{
  return ADVec(
      a[0] - b[0],
      a[1] - b[1],
      a[2] - b[2]);
}

static AD operator*(ADVec a, ADVec b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static AD norm(ADVec a)
{
  return sqrt(a * a);
}

static ADVec cross(ADVec a, ADVec b)
{
  return ADVec(
      a[1] * b[2] - a[2] * b[1],
      a[2] * b[0] - a[0] * b[2],
      a[0] * b[1] - a[1] * b[0]);
}

static AD tet_volume(ADVec v[4])
{
  return (cross((v[1] - v[0]), (v[2] - v[0])) * (v[3] - v[0])) / AD(6);
}

static AD triangle_area(ADVec v[3])
{
  return norm(cross((v[1] - v[0]), (v[2] - v[0]))) / AD(2.);
}

#define PERFECT_TRIANGLE_QUALITY (std::sqrt(3.0) / 4.0)
#define CUBE(x) ((x)*(x)*(x))
#define PERFECT_TET_QUALITY \
  ((std::sqrt(2.0) / 12.0) / CUBE(std::sqrt(PERFECT_TRIANGLE_QUALITY)))

static AD tet_quality(ADVec v[4])
{
  AD sum_asq(0);
  for (int i = 0; i < 4; ++i) {
    ADVec tri_v[3];
    for (int j = 0; j < 3; ++j)
      tri_v[j] = v[apf::tet_tri_verts[i][j]];
    AD a = triangle_area(tri_v);
    sum_asq = sum_asq + (a * a);
  }
  AD arms = sqrt(sum_asq / 4.);
  AD vol = tet_volume(v);
  AD root_arms = sqrt(arms);
  AD quality = vol / CUBE(root_arms);
  return quality / PERFECT_TET_QUALITY;
}

static AD tet_entity_quality(Mesh* m, Entity* tet, Entity* v)
{
  Entity* tv[4];
  m->getDownward(tet, 0, tv);
  ADVec tx[4];
  for (int i = 0; i < 4; ++i) {
    Vector vx;
    m->getPoint(tv[i], 0, vx);
    if (tv[i] == v)
      tx[i] = ADVec(AD(vx[0], 1, 0, 0),
                    AD(vx[1], 0, 1, 0),
                    AD(vx[2], 0, 0, 1));
    else
      tx[i] = ADVec(AD(vx[0]), AD(vx[1]), AD(vx[2]));
  }
  return tet_quality(tx);
}

static AD min_cavity_quality(Mesh* m, apf::Adjacent& tets, Entity* v)
{
  AD min_qual(1.);
  for (size_t i = 0; i < tets.getSize(); ++i) {
    AD tet_qual = tet_entity_quality(m, tets[i], v);
    if (tet_qual.val() < min_qual.val())
      min_qual = tet_qual;
  }
  return min_qual;
}

bool repositionVertex(Mesh* m, Entity* v,
    int max_iters, double initial_speed)
{
  apf::Adjacent tets;
  m->getAdjacent(v, 3, tets);
  double speed = initial_speed;
  AD min_qual = min_cavity_quality(m, tets, v);
  double prev_qual = min_qual.val();
  for (int iter = 0; iter < max_iters; ++iter) {
    prev_qual = min_qual.val();
    double step = (1.0 - min_qual.val()) * speed;
    Vector grad(
        min_qual.dx(0),
        min_qual.dx(1),
        min_qual.dx(2));
    double grad_mag = grad * grad;
    if (!grad_mag)
      break;
    Vector motion = grad * (step / grad_mag);
    Vector vx;
    m->getPoint(v, 0, vx);
    vx += motion;
    m->setPoint(v, 0, vx);
    min_qual = min_cavity_quality(m, tets, v);
    /* diverging, reduce speed */
    if (min_qual.val() < prev_qual)
      speed /= 2;
  }
  return min_qual.val() > 0;
}

}

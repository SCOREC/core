#include <apfMesh.h>
#include <cmath>

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
      a.val() * b.val(),
      (a.dx(0) * b.val() - a.val() * b.dx(0)) / (b.val() * b.val()),
      (a.dx(1) * b.val() - a.val() * b.dx(1)) / (b.val() * b.val()),
      (a.dx(2) * b.val() - a.val() * b.dx(2)) / (b.val() * b.val()));
}

static AD sqrt(AD a)
{
  return AD(
      sqrt(a.val()),
      a.dx(0) / (2. * sqrt(a.val())),
      a.dx(1) / (2. * sqrt(a.val())),
      a.dx(2) / (2. * sqrt(a.val())));
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

#define PERFECT_TRIANGLE_QUALITY (sqrt(3.0) / 4.0)
#define CUBE(x) ((x)*(x)*(x))
#define PERFECT_TET_QUALITY \
  ((sqrt(2.0) / 12.0) / CUBE(sqrt(PERFECT_TRIANGLE_QUALITY)))

AD tet_quality(ADVec v[4])
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

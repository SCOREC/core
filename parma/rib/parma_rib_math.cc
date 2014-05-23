#include "parma_rib_math.h"
#include <complex>

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
    u[1] = std::complex<double>(-1,  sqrt(3)) / 2.;
    u[2] = std::complex<double>(-1, -sqrt(3)) / 2.;
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

static double tr(apf::Matrix3x3 const& m)
{
  return m[0][0] + m[1][1] + m[2][2];
}

/* wikipedia.org/wiki/Characteristic_polynomial#Characteristic_equation */
static void getCharacteristicPolynomial(apf::Matrix3x3 const& A, Cubic& p)
{
  double tA = tr(A);
  double c2 = (1. / 2.) * (tA * tA + tr(A * A));
  p.a = -1;
  p.b = tA;
  p.c = -c2;
  p.d = apf::det(A);
}

static double getPrincipalEigenvalue(apf::Matrix3x3 const& A)
{
  Cubic p;
  getCharacteristicPolynomial(A, p);
  Cubic::Roots l;
  p.solve(l);
  int best = 0;
  for (int i = 1; i < l.n; ++i)
    if (fabs(l.x[i]) > fabs(l.x[best]))
      best = i;
  return l.x[best];
}

static void getEigenvector(apf::Matrix3x3 const& A, double l, apf::Vector3& v)
{
  apf::Matrix3x3 eye(1,0,0,
                     0,1,0,
                     0,0,1);
  apf::Matrix3x3 basis = apf::transpose(A - eye * l);
/* the rows of this ^ matrix should in the general case
   span a plane, or in bad cases a line or nothing.
   We assume the good case and take the most outstanding cross product
   of any pair of rows */
  apf::Vector3 c[3];
  c[0] = apf::cross(basis[0],basis[1]);
  c[1] = apf::cross(basis[1],basis[2]);
  c[2] = apf::cross(basis[2],basis[0]);
  int best = 0;
  for (int i = 1; i < 3; ++i)
    if (c[i].getLength() > c[best].getLength())
      best = i;
/* if the null space is higher than 1-dimensional, this
   will break in div by zero, or divide by something very close to 0 */
  v = c[best].normalize();
}

namespace parma {

void getPrincipalEigenvector(apf::Matrix3x3 const& A, apf::Vector3& v)
{
  double l = getPrincipalEigenvalue(A);
  getEigenvector(A, l, v);
}

}

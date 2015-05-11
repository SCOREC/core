#ifndef APF_GEOMETRY_H
#define APF_GEOMETRY_H

#include "apfVector.h"

namespace apf {

struct Line {
  Vector3 origin;
  Vector3 direction;
  Line(Vector3 const& o, Vector3 const& d);
};

struct Plane {
  Vector3 normal;
  double radius;
  Plane(Vector3 const& n, double r);
  static Plane fromPoints(Vector3 const& a,
      Vector3 const& b, Vector3 const& c);
  double distance(Vector3 const& a) const;
};

bool areClose(double a, double b, double tol);
bool areClose(Vector3 const& a, Vector3 const& b, double tol);
bool areParallel(Vector3 const& a, Vector3 const& b, double tol);
bool areOrthogonal(Vector3 const& a, Vector3 const& b, double tol);

bool areParallel(Plane const& a, Plane const& b, double tol);
bool areClose(Plane const& a, Plane const& b, double tol);
bool areParallel(Line const& a, Plane const& b, double tol);
bool areClose(Line const& a, Plane const& b, double tol);

Line intersect(Plane const& a, Plane const& b);
Vector3 intersect(Line const& a, Plane const& b);

};

#endif

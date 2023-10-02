#ifndef APF_GEOMETRY_H
#define APF_GEOMETRY_H

#include "apfMatrix.h"

namespace apf {

struct Line {
  Vector3 origin;
  Vector3 direction;
  Line();
  Line(Vector3 const& o, Vector3 const& d);
};

struct LineSegment {
  Vector3 start;
  Vector3 end;
  LineSegment();
  LineSegment(Vector3 const& s, Vector3 const& e);
};

struct Plane {
  Vector3 normal;
  double radius;
  Plane(Vector3 const& n, double r);
  static Plane fromPoints(Vector3 const& a,
      Vector3 const& b, Vector3 const& c);
  double distance(Vector3 const& a) const;
};

struct Frame {
  Matrix3x3 linear;
  Vector3 trans;
  Frame();
  Frame(Matrix3x3 const& l, Vector3 const& t);
  static Frame forRotation(Vector3 const& u, double a);
  static Frame forTranslation(Vector3 const& t);
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

Frame operator*(Frame const& a, Frame const& b);
Vector3 operator*(Frame const& a, Vector3 const& b);

double getAngle(Vector3 const& a, Vector3 const& b);

double getDistance(LineSegment const& ls, Vector3 const& p);

/* size is defined as half width(x), half length(y) and half height(z) */
bool withinBox(Vector3 const& point, Vector3 const& center, Vector3 const& size);

bool withinBox(Vector3 const& point, Vector3 const& center, Vector3 const& size,
               Vector3 const& normal1, Vector3 const& normal2);

/* hlen is the half of the length */
bool withinCyl(Vector3 const& point, Vector3 const& center,
               double hlen, double radius);

bool withinCyl(Vector3 const& point, Vector3 const& center,
               double hlen, double radius, Vector3 const& normal);

}

#endif

#include "apfGeometry.h"
#include "apfMatrix.h"
#include <cmath>
#include <cassert>

namespace apf {

Line::Line()
{
}

Line::Line(Vector3 const& o, Vector3 const& d):
  origin(o),
  direction(d.normalize())
{
}

LineSegment::LineSegment()
{
}

LineSegment::LineSegment(Vector3 const& a, Vector3 const& b):
  start(a),
  end(b)
{
}

Plane::Plane(Vector3 const& n, double r):
  normal(n.normalize()),
  radius(r * n.getLength())
{
}

Plane Plane::fromPoints(Vector3 const& a,
      Vector3 const& b, Vector3 const& c)
{
  Vector3 normal = apf::cross(a - c, b - c).normalize();
  double radius = c * normal;
  return Plane(normal, radius);
}

double Plane::distance(Vector3 const& a) const
{
  return normal * a - radius;
}

Frame::Frame()
{
}

Frame::Frame(Matrix3x3 const& l, Vector3 const& t):
  linear(l),
  trans(t)
{
}

Frame Frame::forRotation(Vector3 const& u, double a)
{
  return Frame(rotate(u, a), Vector3(0,0,0));
}

Frame Frame::forTranslation(Vector3 const& t)
{
  return Frame(Matrix3x3(1,0,0,
                         0,1,0,
                         0,0,1),
               t);
}

bool areClose(double a, double b, double tol)
{
  /* there are more advanced and accurate floating
     point comparison algorithms available.
     consider that if this code gets used a lot
     more than it is now.
     "<=" so that using zero tolerance is an exact
     comparison */
  return std::fabs(b - a) <= tol;
}

bool areClose(Vector3 const& a, Vector3 const& b,
    double tol)
{
  return areClose((b - a).getLength(), 0, tol);
}

bool areClose(Plane const& a, Plane const& b, double tol)
{
  /* radius can be zero, hence the parallel check */
  return areParallel(a.normal, b.normal, tol) &&
         areClose(a.normal * a.radius, b.normal * b.radius, tol);
}

bool areClose(Vector3 const& a, Plane const& b, double tol)
{
  return areClose(b.normal * a, b.radius, tol);
}

bool areParallel(Vector3 const& a, Vector3 const& b, double tol)
{
  return areClose(std::fabs(a * b) / (a.getLength() * a.getLength()), 1, tol);
}

bool areOrthogonal(Vector3 const& a, Vector3 const& b, double tol)
{
  return areClose(a * b / (a.getLength() * a.getLength()), 0, tol);
}

bool areParallel(Plane const& a, Plane const& b, double tol)
{
  return areParallel(a.normal, b.normal, tol);
}

bool areParallel(Line const& a, Plane const& b, double tol)
{
  return areOrthogonal(a.direction, b.normal, tol);
}

bool areClose(Line const& a, Plane const& b, double tol)
{
  return areParallel(a, b, tol) && areClose(a.origin, b, tol);
}

Line intersect(Plane const& a, Plane const& b)
{
  /* we have to compute the line
       x = (origin) + (lambda)(direction)
     as an origin,direction pair by combining the first two planes.
     The direction is straightforward: */
  Vector3 direction = cross(a.normal, b.normal);
  /* the origin can be computed by noting that we want
     to satisfy these two equations:
       normal_0 . x = radius_0
       normal_1 . x = radius_1
     and add this arbitrary but stable equation:
       direction . x = 0
     This gives us an AX=B problem. The normals are required
     to be linearly independent unit vectors,
     so A should be invertible */
  Matrix3x3 A;
  A[0] = a.normal;
  A[1] = b.normal;
  A[2] = direction;
  Vector3 B(a.radius, b.radius, 0);
  Vector3 origin = invert(A) * B;
  return Line(origin, direction);
}

Vector3 intersect(Line const& a, Plane const& b)
{
  /* we need to satisfy two equations:
       x = (origin) + (lambda)(direction)     line equation
       (normal).(x) = radius                  plane equation
     this gives us:
       (normal).((origin) + (lambda)(direction)) = radius
       (normal).(origin) + (lambda)((normal).(direction)) = radius
       lambda = (radius - (normal).(origin)) / ((normal).(direction))
     this will only fail if normal.direction=0, which is
     a line parallel to a plane. */
  double denominator = b.normal * a.direction;
  double numerator = b.radius - (b.normal * a.origin);
  double lambda = numerator / denominator;
  return a.origin + (a.direction * lambda);
}

Frame operator*(Frame const& a, Frame const& b)
{
  return Frame(a.linear * b.linear, a.linear * b.trans + a.trans);
}

Vector3 operator*(Frame const& a, Vector3 const& b)
{
  return a.linear * b + a.trans;
}

double getAngle(Vector3 const& a, Vector3 const& b)
{
  return std::acos((a * b) / (a.getLength() * b.getLength()));
}

double getDistance(LineSegment const& ls, Vector3 const& p)
{
  Vector3 direction = (ls.end - ls.start);
  Vector3 ap = p - ls.start;
  Vector3 along = project(ap, direction);
  Vector3 away = reject(ap, direction);
  if (along * direction < 0)
    return ap.getLength();
  if (along.getLength() > direction.getLength())
    return (p - ls.end).getLength();
  return away.getLength();
}

/* assume the box aligned with x-y-z coordinate */
bool withinBox(Vector3 const& point, Vector3 const& center, Vector3 const& size)
{
  assert(size[0]>0 && size[1]>0 && size[2]>0);
  if ( fabs(point[0]-center[0]) <= size[0] &&
       fabs(point[1]-center[1]) <= size[1] &&
	   fabs(point[2]-center[2]) <= size[2] )
	return true;
  else
    return false;
}

/* any alignment of the box */
bool withinBox(Vector3 const& point, Vector3 const& center, Vector3 const& size,
               Vector3 const& normal1, Vector3 const& normal2)
{
  assert(areOrthogonal(normal1, normal2, 0.0));
  Vector3 xAxis = Vector3(1,0,0);
  Vector3 yAxis = Vector3(0,1,0);
  Frame frame1 = Frame::forRotation(normal1, getAngle(normal1,xAxis));
  Frame frame2 = Frame::forRotation(normal2, getAngle(normal2,yAxis));
  return withinBox(frame1*frame2*point, frame1*frame2*center, size);
}

/* assume the cylinder axis along x-axis */
bool withinCyl(Vector3 const& point, Vector3 const& center,
               double hlen, double radius)
{
  assert(hlen>0 && radius>0);
  Vector3 start = center - Vector3(hlen,0,0);
  Vector3 end   = center + Vector3(hlen,0,0);
  if( fabs(point[0]-center[0]) <= hlen &&
      getDistance(LineSegment(start, end),point) <= radius )
	return true;
  else
    return false;
}

/* any alignment of the cylinder */
bool withinCyl(Vector3 const& point, Vector3 const& center,
               double hlen, double radius, Vector3 const& normal)
{
  Vector3 xAxis = Vector3(1,0,0);
  Frame frame = Frame::forRotation(normal, getAngle(normal,xAxis));
  return withinCyl(frame*point, frame*center, hlen, radius);
}

}

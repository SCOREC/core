#include "phBC.h"
#include <apfMatrix.h>
#include <cstdlib>

namespace ph {

struct DebugConstraint
{
  int modelDim;
  int modelTag;
};

}

std::ostream& operator<<(std::ostream& s, ph::DebugConstraint const& dbg)
{
  s << "\nat model entity dim " << dbg.modelDim
    << " tag " << dbg.modelTag << '\n';
  return s;
}

namespace ph {

struct Constraint
{
  Constraint(int n)
  {
    degreesOfFreedom = n;
  }
  virtual ~Constraint() {}
  int degreesOfFreedom;
  virtual void write(int* iBC, double* BC) = 0;
};

static int maxComponent(apf::Vector3 const& v)
{
  double max = 0;
  int max_i = 0;
  for (int i = 0; i < 3; ++i) {
    double mag = fabs(v[i]);
    if (mag > max) {
      max = mag;
      max_i = i;
    }
  }
  return max_i;
}

struct PlaneConstraint : public Constraint
{
  PlaneConstraint():Constraint(2) {}
  apf::Vector3 normal;
  double radius;
  bool operator==(PlaneConstraint const& other) const
  {
    return normal == other.normal && radius == other.radius;
  }
  virtual void write(int* iBC, double* BC)
  {
    int bit = maxComponent(normal) + 3;
    *iBC |= (1<<bit);
    normal.toArray(BC + 3);
    BC[6] = radius;
  }
};

static Constraint* makePlaneConstraint(double* values)
{
  PlaneConstraint* c = new PlaneConstraint();
  c->radius = values[0];
  c->normal.fromArray(values + 1);
/* plane normal must absolutely be a unit vector */
  c->radius *= c->normal.getLength();
  c->normal = c->normal.normalize();
  return c;
}

/* due to the way constraints are written to file,
   the planes forming this constraint are kept around */
struct LineConstraint : public Constraint
{
  LineConstraint():Constraint(1) {}
  PlaneConstraint* planes[2];
  ~LineConstraint()
  {
    delete planes[0];
    delete planes[1];
  }
  virtual void write(int* iBC, double* BC)
  {
    /* okay, we have to set two bits.
       basically, we have to line them up as best
       we can with the two normals.
       there isn't really a good answer to this. */
    /* try max components */
    int comp0 = maxComponent(planes[0]->normal);
    int comp1 = maxComponent(planes[1]->normal);
    if (comp0 == comp1) {
      /* try to resolve a conflict in max components */
      apf::Vector3 a = planes[0]->normal;
      apf::Vector3 b = planes[1]->normal;
      if (b[comp0] > a[comp0])
        std::swap(a,b);
      b[comp0] = 0;
      int comp1 = maxComponent(b);
      assert(comp0 != comp1);
    }
    int bit0 = comp0 + 3;
    int bit1 = comp1 + 3;
    *iBC |= (1<<bit0);
    *iBC |= (1<<bit1);
    planes[0]->normal.toArray(BC + 3);
    BC[6] = planes[0]->radius;
    planes[1]->normal.toArray(BC + 7);
    BC[10] = planes[1]->radius;
  }
};

struct PointConstraint : public Constraint
{
  PointConstraint():Constraint(0) {}
  apf::Vector3 point;
  apf::Vector3 originalDirection_;
  virtual void write(int* iBC, double* BC)
  {
    for (int i = 0; i < 3; ++i) {
      int bit = i + 3;
      *iBC |= (1<<bit);
    }
    double magnitude = point.getLength();
    BC[6] = magnitude;
    if (magnitude) {
      apf::Vector3 direction = point.normalize();
      direction.toArray(BC + 3);
    } else {
      originalDirection_.toArray(BC + 3);
    }
  }
};

static Constraint* makePointConstraint(double* values)
{
  PointConstraint* c = new PointConstraint();
  c->originalDirection_.fromArray(values + 1);
  c->point = c->originalDirection_ * values[0];
  /* just-in-case normalization */
  c->originalDirection_ = c->originalDirection_.normalize();
  return c;
}

static Constraint* takeOne(Constraint* a, Constraint* b, bool takeFirst = true)
{
  if ( ! takeFirst)
    std::swap(a, b);
  delete b;
  return a;
}

/* two point constraints */
static Constraint* combinePoints(Constraint* a, Constraint* b,
    DebugConstraint const& dbg)
{
  PointConstraint* pa = static_cast<PointConstraint*>(a);
  PointConstraint* pb = static_cast<PointConstraint*>(b);
  /* same points, arbitrary victory */
  if (pa->point == pb->point)
    return takeOne(a, b);
  double ma = pa->point.getLength();
  double mb = pb->point.getLength();
  /* any zero magnitude wins (no-slip wins over weaker constraints) */
  if (ma == 0)
    return takeOne(a, b, true);
  if (mb == 0)
    return takeOne(a, b, false);
  /* multiple non-zero point constraints ? we got a problem. */
  std::cerr << "ph error: point overconstrain: ";
  std::cerr << pa->point << " and " << pb->point << dbg;
  abort();
  return 0;
}

static Constraint* combinePlanes(Constraint* a, Constraint* b,
    DebugConstraint const& dbg)
{
  PlaneConstraint* pa = static_cast<PlaneConstraint*>(a);
  PlaneConstraint* pb = static_cast<PlaneConstraint*>(b);
  if (*pa == *pb) /* same plane, arbitrary winner */
    return takeOne(a, b);
  /* the planes are different, so make sure they're not parallel */
  if (pa->normal == pb->normal) {
    std::cerr << "ph error: different parallel planes" << dbg;
    abort();
  }
  /* different intersecting planes, combine into a line constraint */
  LineConstraint* c = new LineConstraint();
  c->planes[0] = pa;
  c->planes[1] = pb;
  return c;
}

static Constraint* combineLinePlane(Constraint* a, Constraint* b,
    DebugConstraint const& dbg)
{
  LineConstraint* pa = static_cast<LineConstraint*>(a);
  PlaneConstraint* pb = static_cast<PlaneConstraint*>(b);
  /* first off, if the plane is the same as one of
     the two that formed the line, we can leave early. */
  for (int i = 0; i < 2; ++i)
    if (*(pa->planes[i]) == *pb)
      return takeOne(a, b, true); /* just keep the line */
  /* okay, we really have 3 distinct planes. lets combine them. */
  /* first, we have to actually compute the line
       x = (origin) + (lambda)(direction)
     as an origin,direction pair by combining the first two planes.
     The direction is straightforward: */
  apf::Vector3 direction = apf::cross(pa->planes[0]->normal,
                                      pa->planes[1]->normal);
  /* the origin can be computed by noting that we want
     to satisfy these two equations:
       normal_0 . x = radius_0
       normal_1 . x = radius_1
     and add this arbitrary but stable equation:
       direction . x = 0
     This gives us an AX=B problem. The normals are required
     to be unit vectors, and combinePlanes ensures they
     are not equal, so A should be invertible */
  apf::Matrix3x3 A;
  A[0] = pa->planes[0]->normal;
  A[1] = pa->planes[1]->normal;
  A[2] = direction;
  apf::Vector3 B(pa->planes[0]->radius,
                 pa->planes[1]->radius,
                 0);
  apf::Vector3 origin = apf::invert(A) * B;
  /* now we just need to satisfy two equations:
       x = (origin) + (lambda)(direction)     line equation
       (normal).(x) = radius                  plane equation
     this gives us:
       (normal).((origin) + (lambda)(direction)) = radius
       (normal).(origin) + (lambda)((normal).(direction)) = radius
       lambda = (radius - (normal).(origin)) / ((normal).(direction))
     this will only fail if normal.direction=0, which is
     a line parallel to a plane. */
  double denominator = pb->normal * direction;
  if (denominator == 0) {
    std::cerr << "ph error: 3 non-intersecting planes" << dbg;
    abort();
  }
  double numerator = pb->radius - (pb->normal * origin);
  double lambda = numerator / denominator;
  PointConstraint* result = new PointConstraint();
  result->point = origin + (direction * lambda);
  result->originalDirection_ = apf::Vector3(1,0,0);
  delete pa;
  delete pb;
  return result;
}

/* combine linear constraints, we assume that
   a is the accumulation and b is a new constraint.
   also accumulation should proceed from most
   to least strict constraint.
   b is expected to be either a point or a plane
   constraint */
static Constraint* combine(Constraint* a, Constraint* b,
    DebugConstraint const& dbg)
{
  if (!a)
    return b;
  assert(a->degreesOfFreedom <= b->degreesOfFreedom);
  if (a->degreesOfFreedom == 0) {
    /* Rule #1: a point constraint beats anything else */
    if (b->degreesOfFreedom != 0)
      return takeOne(a, b, true);
    /* otherwise, compare the two points */
    if (b->degreesOfFreedom == 0)
      return combinePoints(a, b, dbg);
  }
  if (a->degreesOfFreedom == 1 &&
      b->degreesOfFreedom == 2)
    return combineLinePlane(a, b, dbg);
  if (a->degreesOfFreedom == 2 &&
      b->degreesOfFreedom == 2)
    return combinePlanes(a, b, dbg);
  abort();
  return 0;
}

typedef Constraint* (*Make)(double* values);

/* similar flow to getFirstApplied in phBC.cc, but we
   have to account for all first-reachable BCs, not just
   the first one we see. */
Constraint* combineAll(gmi_model* gm, FieldBCs& bcs, Make make,
    gmi_ent* ge, Constraint* a)
{
  double* v = getValuesOn(gm, bcs, ge);
  if (v) {
    DebugConstraint dbg;
    dbg.modelTag = gmi_tag(gm, ge);
    dbg.modelDim = gmi_dim(gm, ge);
    Constraint* b = make(v);
    return combine(a, b, dbg);
  }
  gmi_set* up = gmi_adjacent(gm, ge, gmi_dim(gm, ge) + 1);
  for (int i = 0; i < up->n; ++i)
    a = combineAll(gm, bcs, make, up->e[i], a);
  gmi_free_set(up);
  return a;
}

bool applyVelocityConstaints(gmi_model* gm, BCs& bcs, gmi_ent* e,
    double* BC, int* iBC)
{
  Constraint* c = 0;
  std::string name = "comp3";
  if (hasBC(bcs, name)) {
    FieldBCs& fbcs = bcs.fields[name];
    c = combineAll(gm, fbcs, makePointConstraint, e, c);
  }
  name = "comp1";
  if (hasBC(bcs, name)) {
    FieldBCs& fbcs = bcs.fields[name];
    c = combineAll(gm, fbcs, makePlaneConstraint, e, c);
  }
  if (!c)
    return false;
  c->write(iBC, BC);
  delete c;
  return true;
}

}

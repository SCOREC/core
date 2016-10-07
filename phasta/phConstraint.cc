#include "phBC.h"
#include <apfGeometry.h>
#include <cstdlib>
#include <cassert>
#include <iostream>

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
  PlaneConstraint(apf::Plane const& a):
    Constraint(2),
  plane(a)
  {
  }
  apf::Plane plane;
  bool operator==(PlaneConstraint const& other) const
  {
    return apf::areClose(plane, other.plane, 0.0);
  }
  virtual void write(int* iBC, double* BC)
  {
    int bit = maxComponent(plane.normal) + 3;
    *iBC |= (1<<bit);
    plane.normal.toArray(BC + 3);
    BC[6] = plane.radius;
  }
};

struct PlaneConstraintElas : public Constraint
{
  PlaneConstraintElas(apf::Plane const& a):
    Constraint(2),
  plane(a)
  {
  }
  apf::Plane plane;
  bool operator==(PlaneConstraintElas const& other) const
  {
    return apf::areClose(plane, other.plane, 0.0);
  }
  virtual void write(int* iBC, double* BC)
  {
    int bit = maxComponent(plane.normal) + 14;
    *iBC |= (1<<bit);
    plane.normal.toArray(BC + 16);
    BC[19] = plane.radius;
  }
};

static Constraint* makePlaneConstraint(double* values)
{
  apf::Vector3 normal;
  normal.fromArray(values + 1);
  return new PlaneConstraint(
      apf::Plane(normal, values[0]));
}

static Constraint* makePlaneConstraintElas(double* values)
{
  apf::Vector3 normal;
  normal.fromArray(values + 1);
  return new PlaneConstraintElas(
      apf::Plane(normal, values[0]));
}

/* due to the way constraints are written to file,
   the planes forming this constraint are kept around */
struct LineConstraint : public Constraint
{
  PlaneConstraint* pcs[2];
  LineConstraint():
    Constraint(1)
  {
    pcs[0] = pcs[1] = 0;
  }
  ~LineConstraint()
  {
    delete pcs[0];
    delete pcs[1];
  }
  virtual void write(int* iBC, double* BC)
  {
    /* okay, we have to set two bits.
       basically, we have to line them up as best
       we can with the two normals.
       there isn't really a good answer to this. */
    /* try max components */
    int comp0 = maxComponent(pcs[0]->plane.normal);
    int comp1 = maxComponent(pcs[1]->plane.normal);
    if (comp0 == comp1) {
      /* try to resolve a conflict in max components */
      apf::Vector3 a = pcs[0]->plane.normal;
      apf::Vector3 b = pcs[0]->plane.normal;
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
    pcs[0]->plane.normal.toArray(BC + 3);
    BC[6] = pcs[0]->plane.radius;
    pcs[1]->plane.normal.toArray(BC + 7);
    BC[10] = pcs[1]->plane.radius;
  }
};

struct LineConstraintElas : public Constraint
{
  PlaneConstraintElas* pcs[2];
  LineConstraintElas():
    Constraint(1)
  {
    pcs[0] = pcs[1] = 0;
  }
  ~LineConstraintElas()
  {
    delete pcs[0];
    delete pcs[1];
  }
  virtual void write(int* iBC, double* BC)
  {
    /* okay, we have to set two bits.
       basically, we have to line them up as best
       we can with the two normals.
       there isn't really a good answer to this. */
    /* try max components */
    int comp0 = maxComponent(pcs[0]->plane.normal);
    int comp1 = maxComponent(pcs[1]->plane.normal);
    if (comp0 == comp1) {
      /* try to resolve a conflict in max components */
      apf::Vector3 a = pcs[0]->plane.normal;
      apf::Vector3 b = pcs[0]->plane.normal;
      if (b[comp0] > a[comp0])
        std::swap(a,b);
      b[comp0] = 0;
      int comp1 = maxComponent(b);
      assert(comp0 != comp1);
    }
    int bit0 = comp0 + 14;
    int bit1 = comp1 + 14;
    *iBC |= (1<<bit0);
    *iBC |= (1<<bit1);
    pcs[0]->plane.normal.toArray(BC + 16);
    BC[19] = pcs[0]->plane.radius;
    pcs[1]->plane.normal.toArray(BC + 20);
    BC[23] = pcs[1]->plane.radius;
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

struct PointConstraintElas : public Constraint
{
  PointConstraintElas():Constraint(0) {}
  apf::Vector3 point;
  apf::Vector3 originalDirection_;
  virtual void write(int* iBC, double* BC)
  {
    for (int i = 0; i < 3; ++i) {
      int bit = i + 14;
      *iBC |= (1<<bit);
    }
    double magnitude = point.getLength();
    BC[19] = magnitude;
    if (magnitude) {
      apf::Vector3 direction = point.normalize();
      direction.toArray(BC + 16);
    } else {
      originalDirection_.toArray(BC + 16);
    }
  }
};

struct InterfaceConstraint : public Constraint
{
  
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

static Constraint* makePointConstraintElas(double* values)
{
  PointConstraintElas* c = new PointConstraintElas();
  c->originalDirection_.fromArray(values + 1);
  c->point = c->originalDirection_ * values[0];
  /* just-in-case normalization */
  c->originalDirection_ = c->originalDirection_.normalize();
  return c;
}

static Constraint* takeFirst(Constraint* a, Constraint* b)
{
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
  if (apf::areClose(pa->point, pb->point, 0.0))
    return takeFirst(a, b);
  double ma = pa->point.getLength();
  double mb = pb->point.getLength();
  /* any zero magnitude wins (no-slip wins over weaker constraints) */
  if (ma == 0)
    return takeFirst(a, b);
  if (mb == 0)
    return takeFirst(b, a);
  /* multiple non-zero point constraints ? we got a problem. */
  std::cerr << "ph error: point overconstraint (velocity): ";
  std::cerr << pa->point << " and " << pb->point << dbg;
  abort();
  return 0;
}

static Constraint* combinePointsElas(Constraint* a, Constraint* b,
    DebugConstraint const& dbg)
{
  PointConstraintElas* pa = static_cast<PointConstraintElas*>(a);
  PointConstraintElas* pb = static_cast<PointConstraintElas*>(b);
  /* same points, arbitrary victory */
  if (apf::areClose(pa->point, pb->point, 0.0))
    return takeFirst(a, b);
  double ma = pa->point.getLength();
  double mb = pb->point.getLength();
  /* any zero magnitude wins (no-slip wins over weaker constraints) */
  if (ma == 0)
    return takeFirst(a, b);
  if (mb == 0)
    return takeFirst(b, a);
  /* multiple non-zero point constraints ? we got a problem. */
  std::cerr << "ph error: point overconstraint (mesh-elas): ";
  std::cerr << pa->point << " and " << pb->point << dbg;
  abort();
  return 0;
}

static Constraint* combinePlanes(Constraint* a, Constraint* b,
    DebugConstraint const& dbg)
{
  PlaneConstraint* pa = static_cast<PlaneConstraint*>(a);
  PlaneConstraint* pb = static_cast<PlaneConstraint*>(b);
  if (apf::areClose(pa->plane, pb->plane, 0.0))
    /* same plane, arbitrary winner */
    return takeFirst(a, b);
  /* the planes are different, so make sure they're not parallel */
  if (apf::areParallel(pa->plane, pb->plane, 0.0)) {
    std::cerr << "ph error: different parallel planes (velocity)" << dbg;
    abort();
  }
  /* different intersecting planes, combine into a line constraint */
  LineConstraint* c = new LineConstraint();
  c->pcs[0] = pa;
  c->pcs[1] = pb;
  return c;
}

static Constraint* combinePlanesElas(Constraint* a, Constraint* b,
    DebugConstraint const& dbg)
{
  PlaneConstraintElas* pa = static_cast<PlaneConstraintElas*>(a);
  PlaneConstraintElas* pb = static_cast<PlaneConstraintElas*>(b);
  if (apf::areClose(pa->plane, pb->plane, 0.0))
    /* same plane, arbitrary winner */
    return takeFirst(a, b);
  /* the planes are different, so make sure they're not parallel */
  if (apf::areParallel(pa->plane, pb->plane, 0.0)) {
    std::cerr << "ph error: different parallel planes (mesh-elas)" << dbg;
    abort();
  }
  /* different intersecting planes, combine into a line constraint */
  LineConstraintElas* c = new LineConstraintElas();
  c->pcs[0] = pa;
  c->pcs[1] = pb;
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
    if (apf::areClose(pa->pcs[i]->plane, pb->plane, 0.0))
      return takeFirst(a, b); /* just keep the line */
  /* okay, we really have 3 distinct planes. lets combine them.
     first, we actually compute the line that has been waiting
     to be evaluated until now */
  apf::Line line = apf::intersect(pa->pcs[0]->plane, pa->pcs[1]->plane);
  /* the line may be on the third plane */
  if (apf::areClose(line, pb->plane, 0.0))
    return takeFirst(a, b); /* keep the line */
  /* it may never intersect */
  if (apf::areParallel(line, pb->plane, 0.0)) {
    std::cerr << "line doesn't intersect plane (velocity)" << dbg;
    abort();
  }
  /* okay, there is a legit intersection point. find it. */
  apf::Vector3 point = apf::intersect(line, pb->plane);
  PointConstraint* result = new PointConstraint();
  result->point = point;
  result->originalDirection_ = apf::Vector3(1,0,0);
  delete pa;
  delete pb;
  return result;
}

static Constraint* combineLinePlaneElas(Constraint* a, Constraint* b,
    DebugConstraint const& dbg)
{
  LineConstraintElas* pa = static_cast<LineConstraintElas*>(a);
  PlaneConstraintElas* pb = static_cast<PlaneConstraintElas*>(b);
  /* first off, if the plane is the same as one of
     the two that formed the line, we can leave early. */
  for (int i = 0; i < 2; ++i)
    if (apf::areClose(pa->pcs[i]->plane, pb->plane, 0.0))
      return takeFirst(a, b); /* just keep the line */
  /* okay, we really have 3 distinct planes. lets combine them.
     first, we actually compute the line that has been waiting
     to be evaluated until now */
  apf::Line line = apf::intersect(pa->pcs[0]->plane, pa->pcs[1]->plane);
  /* the line may be on the third plane */
  if (apf::areClose(line, pb->plane, 0.0))
    return takeFirst(a, b); /* keep the line */
  /* it may never intersect */
  if (apf::areParallel(line, pb->plane, 0.0)) {
    std::cerr << "line doesn't intersect plane (mesh-elas)" << dbg;
    abort();
  }
  /* okay, there is a legit intersection point. find it. */
  apf::Vector3 point = apf::intersect(line, pb->plane);
  PointConstraintElas* result = new PointConstraintElas();
  result->point = point;
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
      return takeFirst(a, b);
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

static Constraint* combineElas(Constraint* a, Constraint* b,
    DebugConstraint const& dbg)
{
  if (!a)
    return b;
  assert(a->degreesOfFreedom <= b->degreesOfFreedom);
  if (a->degreesOfFreedom == 0) {
    /* Rule #1: a point constraint beats anything else */
    if (b->degreesOfFreedom != 0)
      return takeFirst(a, b);
    /* otherwise, compare the two points */
    if (b->degreesOfFreedom == 0)
      return combinePointsElas(a, b, dbg);
  }
  if (a->degreesOfFreedom == 1 &&
      b->degreesOfFreedom == 2)
    return combineLinePlaneElas(a, b, dbg);
  if (a->degreesOfFreedom == 2 &&
      b->degreesOfFreedom == 2)
    return combinePlanesElas(a, b, dbg);
  abort();
  return 0;
}

typedef Constraint* (*Make)(double* values);

/* similar flow to getFirstApplied in phBC.cc, but we
   have to account for all first-reachable BCs, not just
   the first one we see. */
Constraint* combineAll(gmi_model* gm, FieldBCs& bcs, Make make,
    gmi_ent* ge, apf::Vector3 const& x, Constraint* a)
{
  double* v = getBCValue(gm, bcs, ge, x);
  if (v) {
    DebugConstraint dbg;
    dbg.modelTag = gmi_tag(gm, ge);
    dbg.modelDim = gmi_dim(gm, ge);
    Constraint* b = make(v);
    return combine(a, b, dbg);
  }
  gmi_set* up = gmi_adjacent(gm, ge, gmi_dim(gm, ge) + 1);
  for (int i = 0; i < up->n; ++i)
    a = combineAll(gm, bcs, make, up->e[i], x, a);
  gmi_free_set(up);
  return a;
}

Constraint* combineAllElas(gmi_model* gm, FieldBCs& bcs, Make make,
    gmi_ent* ge, apf::Vector3 const& x, Constraint* a)
{
  double* v = getBCValue(gm, bcs, ge, x);
  if (v) {
    DebugConstraint dbg;
    dbg.modelTag = gmi_tag(gm, ge);
    dbg.modelDim = gmi_dim(gm, ge);
    Constraint* b = make(v);
    return combineElas(a, b, dbg);
  }
  gmi_set* up = gmi_adjacent(gm, ge, gmi_dim(gm, ge) + 1);
  for (int i = 0; i < up->n; ++i)
    a = combineAllElas(gm, bcs, make, up->e[i], x, a);
  gmi_free_set(up);
  return a;
}

Constraint* combineInterface
(
  gmi_model* gm, FieldBCs& bcs, Make make,
  gmi_ent* ge, apf::Vector3 const& x, Constraint* a
)
{
  double* v = getBCValue(gm, bcs, ge, x);
  if (v) {
    /* The interface attribute only takes an integer value now
       and it is not even used. So, in order to reuse combineElas,
       fake value (u) is used. It is equivalent to mag=0, direction=(1,0,0)
     */
    double u[4] = {0,1,0,0}; 
    DebugConstraint dbg;
    dbg.modelTag = gmi_tag(gm, ge);
    dbg.modelDim = gmi_dim(gm, ge);
    Constraint* b = make(u);
    return combineElas(a, b, dbg);
  }
  gmi_set* up = gmi_adjacent(gm, ge, gmi_dim(gm, ge) + 1);
  for (int i = 0; i < up->n; ++i)
    a = combineInterface(gm, bcs, make, up->e[i], x, a);
  gmi_free_set(up);
  return a;
}

bool applyVelocityConstaints(gmi_model* gm, BCs& bcs, gmi_ent* e,
    apf::Vector3 const& x, double* BC, int* iBC)
{
  Constraint* c = 0;
  std::string name = "comp3";
  if (haveBC(bcs, name)) {
    FieldBCs& fbcs = bcs.fields[name];
    c = combineAll(gm, fbcs, makePointConstraint, e, x, c);
  }
  name = "comp1";
  if (haveBC(bcs, name)) {
    FieldBCs& fbcs = bcs.fields[name];
    c = combineAll(gm, fbcs, makePlaneConstraint, e, x, c);
  }
  if (!c)
    return false;
  c->write(iBC, BC);
  delete c;
  return true;
}

bool applyElasticConstaints(gmi_model* gm, BCs& bcs, gmi_ent* e,
    apf::Vector3 const& x, double* BC, int* iBC)
{
  Constraint* c = 0;
  std::string name;
  name = "DG interface";
  if (haveBC(bcs, name)) {
    FieldBCs& fbcs = bcs.fields[name];
    c = combineInterface(gm, fbcs, makePointConstraintElas, e, x, c);
  }
  name = "comp3_elas";
  if (haveBC(bcs, name)) {
    FieldBCs& fbcs = bcs.fields[name];
    c = combineAllElas(gm, fbcs, makePointConstraintElas, e, x, c);
  }
  name = "comp1_elas";
  if (haveBC(bcs, name)) {
    FieldBCs& fbcs = bcs.fields[name];
    c = combineAllElas(gm, fbcs, makePlaneConstraintElas, e, x, c);
  }
  if (!c)
    return false;
  c->write(iBC, BC);
  delete c;
  return true;
}
}

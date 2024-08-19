/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maSnap.h"
#include "maAdapt.h"
#include "maOperator.h"
#include "maSnapper.h"
#include "maMatchedSnapper.h"
#include "maLayer.h"
#include "maMatch.h"
#include "maDBG.h"
#include <apfGeometry.h>
#include <pcu_util.h>
#include <lionPrint.h>
#include <iostream>
#include <algorithm>

namespace ma {

static size_t isSurfUnderlyingFaceDegenerate(
    apf::Mesh* m,
    Model* g, // this the model entity in question
    int& axis,
    std::vector<double>& values) // sorted in ascending order
{
  int md = m->getModelType(g);
  PCU_ALWAYS_ASSERT(md == 2);

  values.clear();

  Vector bmin;
  Vector bmax;

  m->boundingBox(g, bmin, bmax);
  double tol = 1.0e-6 * (bmax - bmin).getLength();

  bool isPeriodic[2];
  double range[2][2];
  int numPeriodicDims = 0;
  for (int i = 0; i < md; i++) {
    isPeriodic[i] = m->getPeriodicRange(g,i,range[i]);
    if (isPeriodic[i])
      numPeriodicDims++;
    if (range[i][0] > range[i][1])
      std::swap(range[i][0], range[i][1]);
  }

  if (numPeriodicDims != 1) return 0;

  int periodicAxes = isPeriodic[0] ? 0 : 1;
  int degenAxes = 1 - periodicAxes;
  double candidatePeriodicParam =
    (range[periodicAxes][0] + range[periodicAxes][1]) / 2.0;
  for (int i = 0; i < 2; i++) {
    double candidateDegenParam = range[degenAxes][i];
    double param[2];
    Vector uTan;
    Vector vTan;
    param[periodicAxes] = candidatePeriodicParam;
    param[degenAxes] = candidateDegenParam;
    Vector p(param[0], param[1], 0.0);
    m->getFirstDerivative(g, p, uTan, vTan);
    double uTanSize = uTan.getLength();
    double vTanSize = vTan.getLength();
#ifdef HAVE_CAPSTONE
    uTanSize = uTan * uTan;
    vTanSize = vTan * vTan;
#endif
    if (uTanSize < tol || vTanSize < tol) {
      axis = degenAxes;
      values.push_back(candidateDegenParam);
    }
  }

  std::sort(values.begin(), values.end());
  if (values.size() == 2)
    PCU_ALWAYS_ASSERT(values[1] >= values[0]);
  return values.size();
}


/* this is the logic to deal with discontinuous
   periodic parametric coordinates.
   We assume that if the difference between
   the vertex coordinates of a mesh edge is more than half
   the periodic range, then the edge
   crosses over the discontinuity and we need
   to interpolate differently.
   The last parameter "mode" is there since sometimes
   this function needs to be called with mode = 0 and
   other times it needs to be called with mode = 1.
   For example, "interpolateParametricCoordinates".
   Also, this function needs to be called with with
   mode = 0 first. And if the resulting point corresponding
   to the parametric coordinates is not inside the model,
   it should be called again with mode =1 a second time*/
static double interpolateParametricCoordinate(
    double t,
    double a,
    double b,
    double range[2],
    bool isPeriodic,
    int mode)
{
  if ( ! isPeriodic)
    return (1-t)*a + t*b;
  if (range[0] > range[1])
    std::swap(range[0],range[1]);
  if (a > b)
  {
    std::swap(a,b);
    t = 1-t;
  }
  double period = range[1]-range[0];
  double span = b-a;
  if (!mode) {
    if (span < (period/2))
      return (1-t)*a + t*b;
  }
  else {
    if (span <= (period/2))
      return (1-t)*a + t*b;
  }
  a += period;
  double result = (1-t)*a + t*b;
  if (result >= range[1])
    result -= period;
  PCU_ALWAYS_ASSERT(result >= range[0]);
  PCU_ALWAYS_ASSERT(result < range[1]);
  return result;
}

static void interpolateParametricCoordinateOnEdge(
    apf::Mesh* m,
    Model* g,
    double t,
    const Vector& a,
    const Vector& b,
    Vector& p)
{
  double range[2];
  bool isPeriodic = m->getPeriodicRange(g,0,range);
  p[0] = interpolateParametricCoordinate(t,a[0],b[0],range,isPeriodic, 0);
  p[1] = 0.0;
  p[2] = 0.0;

#ifdef HAVE_CAPSTONE
  // account for non-uniform parameterization of model-edge
  Vector X[3];
  Vector para[2] = {a, b};
  for (int i = 0; i < 2; i++)
    m->snapToModel(g, para[i], X[i]);

  double tMax = 1., tMin = 0.;
  m->snapToModel(g, p, X[2]);

  // check if the snap point on the model edge is
  // approximately at same length from either vertices

  double r = (X[0]-X[2]).getLength();
  double s = (X[1]-X[2]).getLength();

  double alpha = t/(1. - t); //parametric ratio

  int num_it = 0;
  while (r/s < 0.95 * alpha || r/s > 1.05 * alpha) {
    if ( r/s > alpha) {
      tMax = t;
      t = (tMin + t) / 2.;
    }
    else {
      tMin = t;
      t = (t + tMax) / 2.;
    }

    p[0] = interpolateParametricCoordinate(t, a[0], b[0], range, isPeriodic, 0);

    m->snapToModel(g, p, X[2]);

    r = (X[0]-X[2]).getLength();
    s = (X[1]-X[2]).getLength();

    if ( num_it > 20) break;
      num_it++;
  }
#endif
}

// convert (phi,theta) on unit sphere to (x,y,z)
// x = sin(theta) cos(phi)
// y = sin(theta) sin(phi)
// z = cos(theta)
// phi in [0 to 2pi]
// theta in [0 to pi]
static Vector getSphere(const Vector& p)
{
  double phi = p[0];
  double the = p[1];
  Vector res(sin(the)*cos(phi),
             sin(the)*sin(phi),
             cos(the));
  return res;
}

// convert (x,y,z) on unit sphere to (theta,phi)
// x = sin(theta) cos(phi)
// y = sin(theta) sin(phi)
// z = cos(theta)
// phi in [0 to 2pi]
// theta in [0 to pi]
static Vector getInvSphere(const Vector& x)
{
  double phi = atan2(x[1], x[0]);
  if (phi < 0)
    phi = phi + 2 * M_PI;
  // the following 3 lines are to avoid errors
  // caused by floating point operations
  double x2 = x[2];
  if (x2 > 1.0) x2 = 1.0;
  if (x2 < -1.0) x2 = -1.0;
  double the = acos(x2);
  Vector res(phi, the, 0.0);
  return res;
}

// map the periodic range [phiMin,phiMax] to [0,2pi]
static double mapPeriodic(double phi, double phiMin, double phiMax){
  return 2. * M_PI * (phi - phiMin) / (phiMax - phiMin);
}

// map the periodic range [0,2pi] back to [phiMin,phiMax]
static double invMapPeriodic(double x, double phiMin, double phiMax){
  return phiMin + x * (phiMax - phiMin) / 2. / M_PI;
}

// map the degenerate range [theMin,theMax] to [0,pi]
// *** for cases where both poles are present.
static double mapDegenerateBoth(double the, double theMin, double theMax)
{
  return M_PI * (the - theMin) / (theMax - theMin);
}

// map the degenerate range [0,pi] back to[theMin,theMax]
// *** for cases where both poles are present.
static double invMapDegenerateBoth(double x, double theMin, double theMax)
{
  return theMin + x * (theMax - theMin) / M_PI;
}

// map the degenerate range [theMin,theMax] to [0,theRef]
// *** for cases where only north  pole is present.
// *** north pole is the pole corresponding to the lower limit
// *** theRef is angle between the following vectors
// *** vectors connecting the center of the degenerate range
// *** vector connecting the center and the pole
// *** vector connecting the center and any point on the trim edge
static double mapDegenerateNorth(double the, double theMin, double theMax, double theRef)
{
  return theRef * (the - theMin) / (theMax - theMin);
}

// map the degenerate range [0,theRef] back to [theMin, theMax]
// *** for cases where only north  pole is present.
// *** north pole is the pole corresponding to the lower limit
static double invMapDegenerateNorth(double x, double theMin, double theMax, double theRef)
{
  return theMin + x * (theMax - theMin) / theRef;
}


// map the degenerate range [theMin,theMax] to [pi-theRef,pi]
// *** for cases where only south  pole is present.
// *** south pole is the pole corresponding to the upper limit limit
// *** theRef is angle between the following vectors
// *** vectors connecting the center of the degenerate range
// *** vector connecting the center and the pole
// *** vector connecting the center and any point on the trim edge
static double mapDegenerateSouth(double the, double theMin, double theMax, double theRef)
{
  return M_PI - theRef + theRef * (the - theMin) / (theMax - theMin);
}

// map the degenerate range [pi-theRef,pi] back to [theMin,theMax]
// *** for cases where only south  pole is present.
// *** south pole is the pole corresponding to the upper limit limit
static double invMapDegenerateSouth(double x, double theMin, double theMax, double theRef)
{
  return theMin + (x + theRef - M_PI) * (theMax - theMin) / theRef;
}

// computes the slerp between to spherical coordinates
// see https://en.wikipedia.org/wiki/Slerp
static Vector getSlerp(double t, Vector a, Vector b)
{
  double omega = acos(a*b/a.getLength()/b.getLength());
  return (a * sin((1.-t)*omega) + b * sin(t*omega)) / sin(omega);
}

static Vector interpolateParametricCoordinatesDoublePoles(
    double t,
    Vector a,
    Vector b,
    double p_range[2],
    double d_range[2],
    int p_axes,
    int d_axes)
{

  double a_the = mapDegenerateBoth(a[d_axes], d_range[0], d_range[1]);
  double a_phi = mapPeriodic(a[p_axes], p_range[0], p_range[1]);
  Vector aTransformed = Vector(a_phi, a_the, 0.0);

  double b_the = mapDegenerateBoth(b[d_axes], d_range[0], d_range[1]);
  double b_phi = mapPeriodic(b[p_axes], p_range[0], p_range[1]);
  Vector bTransformed = Vector(b_phi, b_the, 0.0);

  Vector xa = getSphere(aTransformed);
  Vector xb = getSphere(bTransformed);

  Vector xs = getSlerp(t, xa, xb);

  Vector sTransformed = getInvSphere(xs);

  double s_phi = invMapPeriodic(sTransformed[0], p_range[0], p_range[1]);
  double s_the = invMapDegenerateBoth(sTransformed[1], d_range[0], d_range[1]);

  return (d_axes == 0) ? Vector(s_the, s_phi, 0.0) : Vector(s_phi, s_the, 0.0);
}


static Vector interpolateParametricCoordinatesSinglePole(
    double t,
    Vector a,
    Vector b,
    double p_range[2],
    double d_range[2],
    int p_axes,
    int d_axes,
    double theRef,
    double pp)
{
  double s_phi = 0.0 , s_the = 0.0;
  // check the NORTH pole (i.e., at least one point inside NORTH pole space)
  // definition of NORTH = pp is at the beginning of the d_range
  if (std::abs(pp - d_range[0]) < 1.e-6) {
    double a_the = mapDegenerateNorth(a[d_axes], d_range[0], d_range[1], theRef);
    double a_phi = mapPeriodic(a[p_axes], p_range[0], p_range[1]);
    Vector aTransformed = Vector(a_phi, a_the, 0.0);

    double b_the = mapDegenerateNorth(b[d_axes], d_range[0], d_range[1], theRef);
    double b_phi = mapPeriodic(b[p_axes], p_range[0], p_range[1]);
    Vector bTransformed = Vector(b_phi, b_the, 0.0);

    Vector xa = getSphere(aTransformed);
    Vector xb = getSphere(bTransformed);

    Vector xs = getSlerp(t, xa, xb);

    Vector sTransformed = getInvSphere(xs);

    s_phi = invMapPeriodic(sTransformed[0], p_range[0], p_range[1]);
    s_the = invMapDegenerateNorth(sTransformed[1], d_range[0], d_range[1], theRef);
  }

  // check the SOUTH pole (i.e., at lease one point inside SOUTH pole space)
  // definition of SOUTH = pp is at the end of the d_range
  if (std::abs(pp - d_range[1]) < 1.e-6) {
    double a_the = mapDegenerateSouth(a[d_axes], d_range[0], d_range[1], theRef);
    double a_phi = mapPeriodic(a[p_axes], p_range[0], p_range[1]);
    Vector aTransformed = Vector(a_phi, a_the, 0.0);

    double b_the = mapDegenerateSouth(b[d_axes], d_range[0], d_range[1], theRef);
    double b_phi = mapPeriodic(b[p_axes], p_range[0], p_range[1]);
    Vector bTransformed = Vector(b_phi, b_the, 0.0);

    Vector xa = getSphere(aTransformed);
    Vector xb = getSphere(bTransformed);

    Vector xs = getSlerp(t, xa, xb);

    Vector sTransformed = getInvSphere(xs);

    s_phi = invMapPeriodic(sTransformed[0], p_range[0], p_range[1]);
    s_the = invMapDegenerateSouth(sTransformed[1], d_range[0], d_range[1], theRef);
  }

  return (d_axes == 0) ? Vector(s_the, s_phi, 0.0) : Vector(s_phi, s_the, 0.0);
}

// This handles the following cases:
// 1 - when the surface has both poles
// 2 - when the surface has only one of the poles (i.e., cut sphere)
static void interpolateParametricCoordinatesOnDegenerateFace(
    apf::Mesh* m,
    Model* g,
    double t,
    const Vector& a, // the parametric coords of the 1st point
    const Vector& b, // the parametric coords of the 2nd point
    int d_axes, // the degenerate axes. Can be either 0 or 1
    const std::vector<double>& vals, // singular points on the degenerate axes
    Vector& p)
{
  int numPoles = vals.size();
  PCU_ALWAYS_ASSERT(numPoles == 1 || numPoles == 2);
  /* // convert vectors to uv matrix such that: */
  /* // uv[0][] is the parametric coordinates of the 1st point, and */
  /* // uv[1][] is the parametric coordinates of the 2nd point. */
  /* double uv[2][2]; // = {{a.x(), a.y()} , */
		   /* //   {b.x(), b.y()}}; */
  /* uv[0][0] = a.x(); */
  /* uv[0][1] = b.x(); */
  /* uv[1][0] = a.y(); */
  /* uv[1][1] = b.y(); */
  /* // these are constants */
  /* // TODO: figure out what is the best choice! */
  /* double d_ratio = 0.01; */
  /* double p_ratio = 0.05; */

  double p_range[2];
  double d_range[2];
  int p_axes = 1 - d_axes;

  PCU_ALWAYS_ASSERT(m->getPeriodicRange(g, p_axes, p_range));
  PCU_ALWAYS_ASSERT(!m->getPeriodicRange(g, d_axes, d_range));


  if (numPoles == 2) {
    p = interpolateParametricCoordinatesDoublePoles(t, a, b, p_range, d_range,
    	p_axes, d_axes);
  }

  if (numPoles == 1) {
    Vector x1;
    Vector x2;
    Vector x3;
    Vector p1;
    Vector p2;
    Vector p3;

    if (std::abs(vals[0]-d_range[0]) < 1e-6) { // north pole
      p1 = (d_axes == 0) ? Vector(d_range[0], p_range[0], 0.0) : Vector(p_range[0], d_range[0], 0.0);
      p2 = (d_axes == 0) ? Vector(d_range[1], p_range[0], 0.0) : Vector(p_range[0], d_range[1], 0.0);
      p3 = (d_axes == 0) ? Vector(d_range[1], (p_range[0]+p_range[1])/2, 0.0) :
			   Vector((p_range[0]+p_range[1])/2, d_range[1], 0.0);
    }
    else { // south pole
      p1 = (d_axes == 0) ? Vector(d_range[1], p_range[0], 0.0) : Vector(p_range[0], d_range[1], 0.0);
      p2 = (d_axes == 0) ? Vector(d_range[0], p_range[0], 0.0) : Vector(p_range[0], d_range[0], 0.0);
      p3 = (d_axes == 0) ? Vector(d_range[0], (p_range[0]+p_range[1])/2, 0.0) :
			   Vector((p_range[0]+p_range[1])/2, d_range[0], 0.0);
    }
    m->snapToModel(g, p1, x1);
    m->snapToModel(g, p2, x2);
    m->snapToModel(g, p3, x3);

    double theRef = acos((x2-x1) * (x3-x1) / (x2-x1).getLength() / (x3-x1).getLength());

    p = interpolateParametricCoordinatesSinglePole(t, a, b, p_range, d_range,
    	p_axes, d_axes, theRef, vals[0]);
  }
}

static void interpolateParametricCoordinatesOnRegularFace(
    apf::Mesh* m,
    Model* g,
    double t,
    const Vector& a, // the parametric coords of the 1st point
    const Vector& b, // the parametric coords of the 2nd point
    Vector& p)
{
  double range[2];
  int dim = m->getModelType(g);
  for (int d=0; d < dim; ++d) {
    bool isPeriodic = m->getPeriodicRange(g,d,range);
    p[d] = interpolateParametricCoordinate(t,a[d],b[d],range,isPeriodic, 0);
  }

  /* check if the new point is inside the model.
   * otherwise re-run the above loop with last option
   * in "interpolateParametricCoordinae" being 1.
   * Notes
   * 1) we are assuming manifold surfaces
   * 2) we only check for faces that are periodic
   */

#ifndef HAVE_CAPSTONE
  // this need to be done for faces, only
  if (dim != 2)
    return;

  Vector x;
  bool ok;
  ok = m->isParamPointInsideModel(g, &p[0], x);
  if (ok)
    return;

  for (int d=0; d < dim; ++d) {
    bool isPeriodic = m->getPeriodicRange(g,d,range);
    p[d] = interpolateParametricCoordinate(t,a[d],b[d],range,isPeriodic, 1);
  }
#endif
}

static void interpolateParametricCoordinatesOnFace(
    apf::Mesh* m,
    Model* g,
    double t,
    const Vector& a, // the parametric coords of the 1st point
    const Vector& b, // the parametric coords of the 2nd point
    Vector& p)
{
  std::vector<double> vals;
  int axes;
  size_t num = isSurfUnderlyingFaceDegenerate(m, g, axes, vals);

  if (num > 0) { // the underlying surface is degenerate
#ifndef HAVE_CAPSTONE
    interpolateParametricCoordinatesOnDegenerateFace(m, g, t, a, b, axes, vals, p);
#else
    // account for non-uniform parameterization of model-edge
    Vector X[3];
    Vector para[2] = {a, b};
    for (int i = 0; i < 2; i++)
      m->snapToModel(g, para[i], X[i]);

    interpolateParametricCoordinatesOnDegenerateFace(m, g, t, para[0], para[1], axes, vals, p);

    double tMax = 1., tMin = 0.;
    m->snapToModel(g, p, X[2]);

    // check if the snap point on the model edge is
    // approximately at same length from either vertices

    double r = (X[0]-X[2]).getLength();
    double s = (X[1]-X[2]).getLength();

    double alpha = t/(1. - t); //parametric ratio

    int num_it = 0;
    while (r/s < 0.95 * alpha || r/s > 1.05 * alpha) {
      if ( r/s > alpha) {
	tMax = t;
	t = (tMin + t) / 2.;
      }
      else {
	tMin = t;
	t = (t + tMax) / 2.;
      }

      interpolateParametricCoordinatesOnDegenerateFace(m, g, t, para[0], para[1], axes, vals, p);

      m->snapToModel(g, p, X[2]);

      r = (X[0]-X[2]).getLength();
      s = (X[1]-X[2]).getLength();

      if ( num_it > 20) break;
	num_it++;
    }
#endif
  }
  else
    interpolateParametricCoordinatesOnRegularFace(m, g, t, a, b, p);
}


void interpolateParametricCoordinates(
    apf::Mesh* m,
    Model* g,
    double t,
    Vector const& a,
    Vector const& b,
    Vector& p)
{
  int md = m->getModelType(g);
  if (md == 1)
    interpolateParametricCoordinateOnEdge(m, g, t, a, b, p);
  else if (md == 2)
    interpolateParametricCoordinatesOnFace(m, g, t, a, b, p);
  else
    lion_oprint(1,"model entity must be an edge or a face\n");
}

static void transferParametricBetween(
    Mesh* m,
    Model* g,
    Entity* v[2],
    double t,
    Vector& p)
{
  Vector ep[2];
  for (int i=0; i < 2; ++i)
    m->getParamOn(g,v[i],ep[i]);
  ma::interpolateParametricCoordinates(m,g,t,ep[0],ep[1],p);
}

void transferParametricOnEdgeSplit(
    Mesh* m,
    Entity* e,
    double t,
    Vector& p)
{
  Model* g = m->toModel(e);
  int modelDimension = m->getModelType(g);
  if (m->getDimension()==3 && modelDimension==3) return;
  Entity* ev[2];
  m->getDownward(e,0,ev);
  ma::transferParametricBetween(m, g, ev, t, p);
}

void transferParametricOnTriSplit(
    Mesh* m,
    Entity* face,
    const Vector& xi,
    Vector& param)
{
  Model* g = m->toModel(face);
  int modelDimension = m->getModelType(g);
  if (m->getDimension() == 3 && modelDimension == 3) return;
  Entity* tv[3];
  m->getDownward(face, 0, tv);
  Vector ep[3], pTemp;
  for (int i = 0; i < 3; ++i) {
    m->getParamOn(g, tv[i], ep[i]);
  }
  // TODO: Might be missing some cases here
  for (int d = 0; d < modelDimension; ++d) {
    param[d] = xi[0]*ep[0][d] + xi[1]*ep[1][d] + xi[2]*ep[2][d];
  }
}

void transferParametricOnQuadSplit(
    Mesh* m,
    Entity* quad,
    Entity* v01,
    Entity* v32,
    double y,
    Vector& p)
{
  Model* g = m->toModel(quad);
  int modelDimension = m->getModelType(g);
  if (modelDimension==m->getDimension()) return;
  Entity* v[2];
  v[0] = v01; v[1] = v32;
  ma::transferParametricBetween(m, g, v, y, p);
}

void getClosestPointParametricCoordinates(
    apf::Mesh* m,
    Model* g,
    double t,
    Vector const& a,
    Vector const& b,
    Vector& p)
{
  Vector testPt = a * (1 - t) + b * t;
  Vector targetPt;
  m->getClosestPoint(g, testPt, targetPt, p);
  (void) targetPt;
}

void transferToClosestPointOnEdgeSplit(
    Mesh* m,
    Entity* e,
    double t,
    Vector& p)
{
  Model* g = m->toModel(e);
  int modelDimension = m->getModelType(g);
  if (m->getDimension()==3 && modelDimension==3) return;
  Entity* ev[2];
  m->getDownward(e,0,ev);
  Vector a = getPosition(m, ev[0]);
  Vector b = getPosition(m, ev[1]);
  getClosestPointParametricCoordinates(m, g, t, a, b, p);
}

void transferToClosestPointOnTriSplit(
    Mesh* m,
    Entity* face,
    const Vector& xi,
    Vector& param)
{
  Model* g = m->toModel(face);
  int modelDimension = m->getModelType(g);
  if (m->getDimension() == 3 && modelDimension == 3) return;
  Entity* tv[3];
  m->getDownward(face, 0, tv);
  Vector x[3];

  for (int i = 0; i < 3; i++) {
    x[i] = getPosition(m, tv[i]);
  }

  Vector testPt = x[0] * xi[0] + x[1] * xi[1] + x[2] * xi[2];
  Vector targetPt;

  m->getClosestPoint(g, testPt, targetPt, param);
  (void) targetPt;
}

void transferToClosestPointOnQuadSplit(
    Mesh* m,
    Entity* quad,
    Entity* v01,
    Entity* v32,
    double y,
    Vector& p)
{
  Model* g = m->toModel(quad);
  int modelDimension = m->getModelType(g);
  if (modelDimension==m->getDimension()) return;
  Vector a = getPosition(m, v01);
  Vector b = getPosition(m, v32);
  Vector testPt = a * (1 - y) + b * y;
  Vector targetPt;
  m->getClosestPoint(g, testPt, targetPt, p);
  (void) targetPt;
}

static void getSnapPoint(Mesh* m, Entity* v, Vector& x)
{
  m->getPoint(v,0,x);
  Vector p;
  m->getParam(v,p);
  Model* g = m->toModel(v);
  m->snapToModel(g,p,x);
}

class SnapAll : public Operator
{
  public:
    SnapAll(Adapt* a, Tag* t, bool simple):
      snapper(a, t, simple)
    {
      adapter = a;
      tag = t;
      successCount = 0;
      didAnything = false;
      vert = 0;
    }
    int getTargetDimension() {return 0;}
    bool shouldApply(Entity* e)
    {
      if ( ! getFlag(adapter, e, SNAP))
        return false;
      vert = e;
      snapper.setVert(e);
      return true;
    }
    bool requestLocality(apf::CavityOp* o)
    {
      return snapper.requestLocality(o);
    }
    void apply()
    {
      bool snapped = snapper.run();
      didAnything = didAnything || snapped || snapper.dug;
      if (snapped)
        ++successCount;
      clearFlag(adapter, vert, SNAP);
    }
    int successCount;
    bool didAnything;
  private:
    Adapt* adapter;
    Tag* tag;
    Entity* vert;
    Snapper snapper;
};

bool snapAllVerts(Adapt* a, Tag* t, bool isSimple, long& successCount)
{
  SnapAll op(a, t, isSimple);
  applyOperator(a, &op);
  successCount += a->mesh->getPCU()->Add<long>(op.successCount);
  return a->mesh->getPCU()->Or(op.didAnything);
}

class SnapMatched : public Operator
{
  public:
    SnapMatched(Adapt* a, Tag* t, bool simple):
      snapper(a, t, simple)
    {
      adapter = a;
      tag = t;
      successCount = 0;
      didAnything = false;
      vert = 0;
    }
    int getTargetDimension() {return 0;}
    bool shouldApply(Entity* e)
    {
      if ( ! getFlag(adapter, e, SNAP))
        return false;
      vert = e;
      snapper.setVert(e);
      return true;
    }
    bool requestLocality(apf::CavityOp* o)
    {
      return snapper.requestLocality(o);
    }
    void apply()
    {
      snapper.setVerts();
      bool snapped = snapper.trySnaps();
      didAnything = didAnything || snapped;
      if (snapped)
        ++successCount;
      clearFlagMatched(adapter, vert, SNAP);
    }
    int successCount;
    bool didAnything;
  private:
    Adapt* adapter;
    Tag* tag;
    Entity* vert;
    MatchedSnapper snapper;
};

bool snapMatchedVerts(Adapt* a, Tag* t, bool isSimple, long& successCount)
{
  SnapMatched op(a, t, isSimple);
  applyOperator(a, &op);
  successCount += a->mesh->getPCU()->Add<long>(op.successCount);
  return a->mesh->getPCU()->Or(op.didAnything);
}

long tagVertsToSnap(Adapt* a, Tag*& t)
{
  Mesh* m = a->mesh;
  int dim = m->getDimension();
  t = m->createDoubleTag("ma_snap", 3);
  Entity* v;
  long n = 0;
  Iterator* it = m->begin(0);
  while ((v = m->iterate(it))) {
    int md = m->getModelType(m->toModel(v));
    if (dim == 3 && md == 3)
      continue;
    Vector s;
    getSnapPoint(m, v, s);
    Vector x = getPosition(m, v);
    if (apf::areClose(s, x, 1e-12))
      continue;
    m->setDoubleTag(v, t, &s[0]);
    if (m->isOwned(v))
      ++n;
  }
  m->end(it);
  return m->getPCU()->Add<long>(n);
}

static void markVertsToSnap(Adapt* a, Tag* t)
{
  HasTag p(a->mesh, t);
  markEntities(a, 0, p, SNAP, DONT_SNAP);
}

bool snapOneRound(Adapt* a, Tag* t, bool isSimple, long& successCount)
{
  markVertsToSnap(a, t);
  if (a->mesh->hasMatching())
    return snapMatchedVerts(a, t, isSimple, successCount);
  else
    return snapAllVerts(a, t, isSimple, successCount);
}

long snapTaggedVerts(Adapt* a, Tag* tag)
{
  long successCount = 0;
  /* there are two approaches possible here:
   * 1- first snap all the vertices we can without any additional
   * operation such as digging (simple snap). And then try snapping
   * the remaining vertices that will need extra modifications (non-
   * simple snap).
   * 2- first do the non-simple snaps and then the simple snaps.
   *
   * Here we choose approach 2 for the following reasons
   * (a) approach 2 is approximately as fast as approach 1
   * (b) the problematic snaps will be attempted as soon as possible.
   * This is extremely helpful because if we wait until later on
   * bringing the vert to-be-snapped to the boundary might become more
   * difficult due to the change in location of neighboring verticies
   * that will be snapped before the problematic vert to-be-snapped.
   */
  while (snapOneRound(a, tag, false, successCount));
  while (snapOneRound(a, tag, true, successCount));
  return successCount;
}

void snap(Adapt* a)
{
  if ( ! a->input->shouldSnap)
    return;
  double t0 = pcu::Time();
  Tag* tag;
  /* we are starting to support a few operations on matched
     meshes, including snapping+UR. this should prevent snapping
     from modifying any matched entities */
  preventMatchedCavityMods(a);
  long targets = tagVertsToSnap(a, tag);
  long success = snapTaggedVerts(a, tag);
  snapLayer(a, tag);
  apf::removeTagFromDimension(a->mesh, tag, 0);
  a->mesh->destroyTag(tag);
  double t1 = pcu::Time();
  print(a->mesh->getPCU(), "snapped in %f seconds: %ld targets, %ld non-layer snaps",
    t1 - t0, targets, success);
  if (a->hasLayer)
    checkLayerShape(a->mesh, "after snapping");
}

void visualizeGeometricInfo(Mesh* m, const char* name)
{
  Tag* dimensionTag = m->createIntTag("ma_geom_dim",1);
  Tag* idTag = m->createIntTag("ma_geom_id",1);
  apf::Field* field = apf::createLagrangeField(m,"ma_param",apf::VECTOR,1);
  apf::Field* targetField = apf::createLagrangeField(m,"ma_target",apf::VECTOR,1);
  Iterator* it = m->begin(0);
  Entity* v;
  while ((v = m->iterate(it)))
  {
    Model* c = m->toModel(v);
    int dimension = m->getModelType(c);
    m->setIntTag(v,dimensionTag,&dimension);
    int id = m->getModelTag(c);
    m->setIntTag(v,idTag,&id);
    Vector p;
    Vector xp = getPosition(m, v);
    m->getParam(v,p);
    if (dimension == 2 || dimension == 1) {
      Vector x;
      m->isParamPointInsideModel(c, &p[0], x);
      apf::setVector(targetField, v, 0, x - xp);
    }
    else {
      Vector x(0., 0., 0.);
      apf::setVector(targetField, v, 0, x);
    }
    apf::setVector(field,v,0,p);
  }
  m->end(it);
  apf::writeVtkFiles(name,m);
  it = m->begin(0);
  while ((v = m->iterate(it)))
  {
    m->removeTag(v,dimensionTag);
    m->removeTag(v,idTag);
  }
  m->end(it);
  m->destroyTag(dimensionTag);
  m->destroyTag(idTag);
  apf::destroyField(field);
  apf::destroyField(targetField);
}

}

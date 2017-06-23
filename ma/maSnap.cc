/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include <PCU.h>
#include "maSnap.h"
#include "maAdapt.h"
#include "maOperator.h"
#include "maSnapper.h"
#include "maLayer.h"
#include "maMatch.h"
#include <apfGeometry.h>
#include <pcu_util.h>
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


  double tol = 1.0e-12;
  values.clear();

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
    m->getFirstDerivative(g, param, uTan, vTan);
    double uTanSize = uTan.getLength();
    double vTanSize = vTan.getLength();
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
}


static Vector interpolateParametricCoordinatesDoublePoles(
    double t,
    double uv[2][2],
    double p_range[2],
    int p_axes,
    int d_axes,
    double p_ratio,
    double d_ratio,
    double sp,
    double np)
{
  double d_threshold, p_threshold, p_size, d_size;
  p_size = std::abs(p_range[1] - p_range[0]);
  d_size = std::abs(np - sp);
  p_threshold = p_ratio * p_size;
  d_threshold = d_ratio * d_size;

  // From now on
  // the(ta) is the degenerate coordinate
  // phi     is the periodic   coordinate

  // the corresponding interpolated values
  double the_t = 0.0;
  double phi_t = 0.0;
  double phi_1 = uv[p_axes][0];
  double phi_2 = uv[p_axes][1];
  double the_1 = uv[d_axes][0];
  double the_2 = uv[d_axes][1];

  double phi_span = std::abs(phi_2 - phi_1);

  // First, take care of the cases where one or both of the
  // thetas (degenerate coord) are close to the poles
  bool flag = false;
  // sort points such that the_1 < the_2
  if (the_1 > the_2) {
    std::swap(phi_1, phi_2);
    std::swap(the_1, the_2);
    t = 1 - t;
  }
  the_t = (1 - t) * the_1 + t * the_2;
  double phi_tmp = interpolateParametricCoordinate(t, phi_1, phi_2, p_range, true, 0);

  // check the south pole (i.e., at least one point inside south pole space)
  if (std::abs(the_1 - sp) <= d_threshold) {
    if (std::abs(the_2 - sp) <= d_threshold) {
      if (std::abs(phi_span - p_size/2.) > p_threshold) {
	phi_t = phi_tmp;
	flag = true;
      }
    }
    else {
      phi_t = phi_2;
      flag = true;
    }
  }
  // check the north pole (i.e., at lease one point inside north pole space)
  if (std::abs(the_2 - np) <= d_threshold) {
    if (std::abs(the_1 - np) <= d_threshold) {
      if (std::abs(phi_span - p_size/2.) > p_threshold) {
	phi_t = phi_tmp;
	flag = true;
      }
    }
    else {
      phi_t = phi_1;
      flag = true;
    }
  }
  // check the case when neither of the points falls in any of the poles
  if (std::abs(the_2 - np) > d_threshold &&
      std::abs(the_1 - sp) > d_threshold) {
    if (std::abs(phi_span - p_size/2.) > p_threshold) {
      phi_t = phi_tmp;
      flag = true;
    }
  }


  // Second, take care of the cases where the connecting line should
  // pass over a pole.
  if (!flag && std::abs(phi_span - p_size/2.) <= p_threshold) {
    // sort points such that phi_1 < phi_2
    if (phi_1 > phi_2) {
      std::swap(phi_1, phi_2);
      std::swap(the_1, the_2);
      t = 1 - t;
    }

    double the_1p, the_2p, the_tp, pole;
    if (the_1 >= 0) {
      pole = np;
      if (the_1 + the_2 >= 0) {
	the_1p = the_1;
	the_2p = the_2;
	the_tp = (1 - t) * the_1p + t * (d_size - the_2p);
	the_t = (the_tp <= pole) ? the_tp : d_size - the_tp;
	phi_t = (the_tp <= pole) ? phi_1 : phi_2;
      }
      else { // the_1 + the_2 < 0
	the_1p =  the_1;
	the_2p = -the_2;
	the_tp = (1 - t) * the_1p + t * (d_size - the_2p);
	the_t = (the_tp <= pole) ? the_tp : d_size - the_tp;
	the_t = -the_t;
	phi_t = (the_tp <= pole) ? phi_1 : phi_2;
      }
    }
    else { // the_1 < 0
      pole = sp;
      if (the_1 + the_2 <= 0) {
	the_1p = the_1;
	the_2p = the_2;
	the_tp = (1 - t) * the_1p + t * (-d_size - the_2p);
	the_t = (the_tp >= pole) ? the_tp : -d_size - the_tp;
	phi_t = (the_tp >= pole) ? phi_1 : phi_2;
      }
      else { // the_1 + the_2 > 0
	the_1p =  the_1;
	the_2p = -the_2;
	the_tp = (1 - t) * the_1p + t * (-d_size - the_2p);
	the_t = (the_tp >= pole) ? the_tp : -d_size - the_tp;
	the_t = -the_t;
	phi_t = (the_tp >= pole) ? phi_1 : phi_2;
      }
    }
  }
  return (d_axes == 0) ? Vector(the_t, phi_t, 0.0) : Vector(phi_t, the_t, 0.0);
}


static Vector interpolateParametricCoordinatesSinglePole(
    double t,
    double uv[2][2],
    double p_range[2],
    int p_axes,
    int d_axes,
    double p_ratio,
    double d_ratio,
    double pp)
{
  double d_threshold, p_threshold, p_size, d_size;
  p_size = std::abs(p_range[1] - p_range[0]);
  d_size = std::abs(2.0 * pp);
  p_threshold = p_ratio * p_size;
  d_threshold = d_ratio * d_size;

  // From now on
  // the(ta) is the degenerate coordinate
  // phi     is the periodic   coordinate

  // the corresponding interpolated values
  double the_t = 0.0;
  double phi_t = 0.0;
  double phi_1 = uv[p_axes][0];
  double phi_2 = uv[p_axes][1];
  double the_1 = uv[d_axes][0];
  double the_2 = uv[d_axes][1];

  double phi_span = std::abs(phi_2 - phi_1);

  // First, take care of the cases where one or both of the
  // thetas (degenerate coord) are close to the poles
  bool flag = false;
  // sort points such that the_1 < the_2
  if (the_1 > the_2) {
    std::swap(phi_1, phi_2);
    std::swap(the_1, the_2);
    t = 1 - t;
  }
  the_t = (1 - t) * the_1 + t * the_2;
  double phi_tmp = interpolateParametricCoordinate(t, phi_1, phi_2, p_range, true, 0);

  // check the south pole (i.e., at least one point inside south pole space)
  if (pp < 0) {
    double sp = pp;
    if (std::abs(the_1 - sp) <= d_threshold) {
      if (std::abs(the_2 - sp) <= d_threshold) {
	if (std::abs(phi_span - p_size/2.) > p_threshold) {
	  phi_t = phi_tmp;
	  flag = true;
	}
      }
      else {
	phi_t = phi_2;
	flag = true;
      }
    }
    else {
      if (std::abs(phi_span - p_size/2.) > p_threshold) {
	phi_t = phi_tmp;
	flag = true;
      }
    }
  }
  // check the north pole (i.e., at lease one point inside north pole space)
  if (pp > 0) {
    double np = pp;
    if (std::abs(the_2 - np) <= d_threshold) {
      if (std::abs(the_1 - np) <= d_threshold) {
	if (std::abs(phi_span - p_size/2.) > p_threshold) {
	  phi_t = phi_tmp;
	  flag = true;
	}
      }
      else {
	phi_t = phi_1;
	flag = true;
      }
    }
    else {
      if (std::abs(phi_span - p_size/2.) > p_threshold) {
	phi_t = phi_tmp;
	flag = true;
      }
    }
  }

  // Second, take care of the cases where the connecting line should
  // pass over a pole.
  if (!flag && std::abs(phi_span - p_size/2.) <= p_threshold) {
    // sort points such that phi_1 < phi_2
    if (phi_1 > phi_2) {
      std::swap(phi_1, phi_2);
      std::swap(the_1, the_2);
      t = 1 - t;
    }

    double the_1p, the_2p, the_tp, pole;
    if (pp > 0) {
      pole = pp;
      if (the_1 + the_2 >= 0) {
	the_1p = the_1;
	the_2p = the_2;
	the_tp = (1 - t) * the_1p + t * (d_size - the_2p);
	the_t = (the_tp <= pole) ? the_tp : d_size - the_tp;
	phi_t = (the_tp <= pole) ? phi_1 : phi_2;
      }
      else { // the_1 + the_2 < 0
	the_1p =  the_1;
	the_2p = -the_2;
	the_tp = (1 - t) * the_1p + t * (d_size - the_2p);
	the_t = (the_tp <= pole) ? the_tp : d_size - the_tp;
	the_t = -the_t;
	phi_t = (the_tp <= pole) ? phi_1 : phi_2;
      }
    }
    if (pp < 0) {
      pole = pp;
      if (the_1 + the_2 <= 0) {
	the_1p = the_1;
	the_2p = the_2;
	the_tp = (1 - t) * the_1p + t * (-d_size - the_2p);
	the_t = (the_tp >= pole) ? the_tp : -d_size - the_tp;
	phi_t = (the_tp >= pole) ? phi_1 : phi_2;
      }
      else { // the_1 + the_2 > 0
	the_1p =  the_1;
	the_2p = -the_2;
	the_tp = (1 - t) * the_1p + t * (-d_size - the_2p);
	the_t = (the_tp >= pole) ? the_tp : -d_size - the_tp;
	the_t = -the_t;
	phi_t = (the_tp >= pole) ? phi_1 : phi_2;
      }
    }
  }
  return (d_axes == 0) ? Vector(the_t, phi_t, 0.0) : Vector(phi_t, the_t, 0.0);
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
  // convert vectors to uv matrix such that:
  // uv[0][] is the parametric coordinates of the 1st point, and
  // uv[1][] is the parametric coordinates of the 2nd point.
  double uv[2][2]; // = {{a.x(), a.y()} ,
		   //   {b.x(), b.y()}};
  uv[0][0] = a.x();
  uv[0][1] = b.x();
  uv[1][0] = a.y();
  uv[1][1] = b.y();
  // these are constants
  // TODO: figure out what is the best choice!
  double d_ratio = 0.01;
  double p_ratio = 0.05;

  double p_range[2];
  int p_axes = 1 - d_axes;

  PCU_ALWAYS_ASSERT(m->getPeriodicRange(g, p_axes, p_range));

  if (numPoles == 2) {
    p = interpolateParametricCoordinatesDoublePoles(t, uv, p_range,
    	p_axes, d_axes, p_ratio, d_ratio, vals[0], vals[1]);
  }

  if (numPoles == 1) {
    p = interpolateParametricCoordinatesSinglePole(t, uv, p_range,
    	p_axes, d_axes, p_ratio, d_ratio, vals[0]);
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
  if (num > 0) // the underlying surface is degenerate
    interpolateParametricCoordinatesOnDegenerateFace(m, g, t, a, b, axes, vals, p);
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
    printf("model entity must be an edge or a face\n");
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
  if (modelDimension==m->getDimension()) return;
  Entity* ev[2];
  m->getDownward(e,0,ev);
  ma::transferParametricBetween(m, g, ev, t, p);
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
      return true;
    }
    bool requestLocality(apf::CavityOp* o)
    {
      return snapper.setVert(vert, o);
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
    if (md == dim)
      continue;
    Vector s;
    getSnapPoint(m, v, s);
    Vector x = getPosition(m, v);
    if (apf::areClose(s, x, 0.0))
      continue;
    m->setDoubleTag(v, t, &s[0]);
    if (m->isOwned(v))
      ++n;
  }
  m->end(it);
  return PCU_Add_Long(n);
}

static void markVertsToSnap(Adapt* a, Tag* t)
{
  HasTag p(a->mesh, t);
  markEntities(a, 0, p, SNAP, DONT_SNAP);
}

bool snapOneRound(Adapt* a, Tag* t, bool isSimple, long& successCount)
{
  markVertsToSnap(a, t);
  SnapAll op(a, t, isSimple);
  applyOperator(a, &op);
  successCount += PCU_Add_Long(op.successCount);
  return PCU_Or(op.didAnything);
}

long snapTaggedVerts(Adapt* a, Tag* tag)
{
  long successCount = 0;
  /* first snap all the vertices we can without digging.
     This is fast because it uses just the elements around
     the vertex and doesn't think much, it should also handle
     the vast majority of vertices */
  while (snapOneRound(a, tag, true, successCount));
  /* all the remaining vertices now need some kind of modification
     in order to snap.
     Here we turn on the "try digging before snapping" flag,
     which requires two-layer cavities so hopefully fewer vertices
     are involved here */
  while (snapOneRound(a, tag, false, successCount));
  return successCount;
}

void snap(Adapt* a)
{
  if ( ! a->input->shouldSnap)
    return;
  double t0 = PCU_Time();
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
  double t1 = PCU_Time();
  print("snapped in %f seconds: %ld targets, %ld non-layer snaps",
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

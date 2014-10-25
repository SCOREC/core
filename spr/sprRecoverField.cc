/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "spr.h"

#include <PCU.h>

#include <apfMesh.h>
#include <apfShape.h>
#include <apfCavityOp.h>

#include <set>

namespace spr {

/* overall information useful during recovery */
struct Recovery {
  apf::Mesh* mesh;
  /* order of output field and fit polynomial.
     so far limited to 1 or 2,
     and determined by the order of the mesh's
     coordinate field. */
  int order;
  /* the number of terms in a polynomial of the above order
     defined over 3D space: f(x,y,z) */
  int polynomial_terms;
  /* the number of integration points per element.
     this immediately assumes the mesh has one element type.
     warning: the number of points per element is not necessarily
     determined by the order of the output field !
     users sometimes give us "too much" information, i.e. 5 integration
     points per tet to derive a linear field. this is ok, the
     way this is programmed it should handle any nonzero number of points
     per element regardless of output order (by growing bigger patches) */
  int points_per_element;
  /* input field containing integration point data for all elements */
  apf::Field* f;
  /* output field containing recovered nodal data */
  apf::Field* f_star;
};

static int determinePointsPerElement(apf::Field* f)
{
  apf::Mesh* m = apf::getMesh(f);
  int element_type = apf::getFirstType(m, m->getDimension());
  apf::FieldShape* s = apf::getShape(f);
  return s->countNodesOn(element_type);
}

static apf::Field* makeRecoveredField(Recovery* r)
{
  std::string name = "spr_";
  name += apf::getName(r->f);
  return apf::createLagrangeField(
        r->mesh, name.c_str(), apf::getValueType(r->f), r->order);
}

static int countPolynomialTerms(int order)
{
  return ((order + 1) * (order + 2) * (order + 3)) / 6;
}

static void setupRecovery(Recovery* r, apf::Field* f)
{
  r->mesh = apf::getMesh(f);
  r->order = r->mesh->getShape()->getOrder();
  r->polynomial_terms = countPolynomialTerms(r->order);
  r->points_per_element = determinePointsPerElement(f);
  r->f = f;
  r->f_star = makeRecoveredField(r);
}

typedef std::set<apf::MeshEntity*> EntitySet;

struct Patch {
  apf::Mesh* mesh;
  Recovery* recovery;
  /* the entity around which the patch
     is centered. a patch collects elements
     around this entity and then their integration
     points will be used to recover values for
     all nodes on this entity */
  apf::MeshEntity* entity;
  EntitySet elements;
};

static void setupPatch(Patch* p, Recovery* r)
{
  p->mesh = r->mesh;
  p->recovery = r;
  p->entity = 0;
}

static void startPatch(Patch* p, apf::MeshEntity* e)
{
  p->elements.clear();
  p->entity = e;
}

static int countPatchPoints(Patch* p)
{
  return p->recovery->points_per_element * p->elements.size();
}

static bool hasEnoughPoints(Patch* p)
{
  return countPatchPoints(p) >= p->recovery->polynomial_terms;
}

void addElementToPatch(Patch* p, apf::MeshEntity* e)
{
  p->elements.insert(e);
}

void addElementsToPatch(Patch* p, apf::DynamicArray<apf::MeshEntity*>& es)
{
  for (std::size_t i=0; i < es.getSize(); ++i)
    addElementToPatch(p, es[i]);
}

static bool getInitialPatch(Patch* p, apf::CavityOp* o)
{
  if ( ! o->requestLocality(&p->entity,1))
    return false;
  apf::DynamicArray<apf::MeshEntity*> adjacent;
  p->mesh->getAdjacent(p->entity, p->mesh->getDimension(), adjacent);
  addElementsToPatch(p, adjacent);
  return true;
}

static bool addElementsThatShare(Patch* p, int dim,
    EntitySet& old_elements, apf::CavityOp* o)
{
  EntitySet bridges;
  APF_ITERATE(EntitySet, old_elements, it)
  {
    apf::Downward down;
    int nd = p->mesh->getDownward(*it, dim, down);
    for (int i=0; i < nd; ++i)
      bridges.insert(down[i]);
  }
  std::vector<apf::MeshEntity*> 
    bridge_array(bridges.begin(),bridges.end());
  bridges.clear();
  if ( ! o->requestLocality(&(bridge_array[0]),bridge_array.size()))
    return false;
  for (size_t i=0; i < bridge_array.size(); ++i)
  {
    apf::Adjacent candidates;
    p->mesh->getAdjacent(bridge_array[i], p->mesh->getDimension(), candidates);
    addElementsToPatch(p, candidates);
  }
  return true;
}

static bool expandAsNecessary(Patch* p, apf::CavityOp* o)
{
  if (hasEnoughPoints(p))
    return true;
  EntitySet old_set = p->elements;
  int d = p->mesh->getDimension();
  for (int shared_dim = d-1; shared_dim >= 0; --shared_dim)
  {
    if (!addElementsThatShare(p, shared_dim, old_set, o))
      return false;
    if (hasEnoughPoints(p));
      return true;
  }
  bool hope = p->elements.size() > old_set.size();
  if (hope)
    return expandAsNecessary(p, o);
  else
  {
    apf::fail("SPR: patch construction: all hope is lost.");
    return false;
  }
}

static bool buildPatch(Patch* p, apf::CavityOp* o)
{
  if (!getInitialPatch(p, o)) return false;
  if (!expandAsNecessary(p, o)) return false;
  return true;
}

struct Samples {
  Samples():num_points(0) {}
  void allocate(int np, int nc)
  {
    num_points = np;
    points.allocate(np);
    values.allocate(np);
    for (int i=0; i < np; ++i)
      values[i].allocate(nc);
  }
  int num_points;
  apf::NewArray<apf::Vector3> points;
  apf::NewArray<apf::NewArray<double> > values;
};

/** @brief get spr point data from element patch
  * @details assumes constant #IP/element
  */
void getSamples(Patch* p, Samples* s)
{
  Recovery* r = p->recovery;
  int np = countPatchPoints(p);
  int nc = apf::countComponents(r->f);
  s->allocate(np,nc);
  std::size_t i = 0;
  APF_ITERATE(EntitySet, p->elements, it) {
    apf::MeshElement* me = apf::createMeshElement(r->mesh, *it);
    for (int l = 0; l < r->points_per_element; ++l) {
      apf::Vector3 param;
      apf::getIntPoint(me, r->order, l, param);
      apf::mapLocalToGlobal(me, param, s->points[i]);
      apf::getComponents(r->f, *it, l, &(s->values[i][0]));
      ++i;
    }
    apf::destroyMeshElement(me);
  }
}

void evalPolynomialTerms(int order,
                         apf::Vector3 const& point,
                         apf::DynamicVector& terms)
{
  apf::Vector3 const& x = point;
  if (order == 1)
  {
    terms.setSize(4);
    terms(0) = 1.0;
    terms(1) = x[0];
    terms(2) = x[1];
    terms(3) = x[2];
  }
  else if (order == 2)
  {
    terms.setSize(10);
    terms(0) = 1.0;
    terms(1) = x[0];
    terms(2) = x[1];
    terms(3) = x[2];
    terms(4) = x[0]*x[1];
    terms(5) = x[1]*x[2];
    terms(6) = x[2]*x[0];
    terms(7) = x[0]*x[0];
    terms(8) = x[1]*x[1];
    terms(9) = x[2]*x[2];
  }
  else
    apf::fail("SPR: invalid polynomial order");
}

/** @brief least squares solution of polynomial coefficients
  * @param order (In) order of polynomial to be fit to data
  * @param num_points (In) number of integration points in patch
  * @param points (In) coordinates of integration points in patch
  * @param values (In) data value associated with integration points
  * @param coeffs (Out) coefficients of fit polynomial
  */
void fitPolynomial(int order,
                   int num_points,
                   apf::NewArray<apf::Vector3> const& points,
                   apf::NewArray<double> const& values,
                   apf::DynamicVector& coeffs)
{
  int m = num_points;
  int n = countPolynomialTerms(order);
  assert(m >= n);
  apf::DynamicMatrix P(m,n);
  apf::DynamicVector p;
  for (int i = 0; i < m; ++i)
  {
    evalPolynomialTerms(order,points[i],p);
    P.setRow(i,p);
  }
  apf::DynamicMatrix PT;
  apf::transpose(P,PT);
  apf::DynamicMatrix A;
  apf::multiply(PT,P,A);
  apf::DynamicVector v(m);
  for (int i=0; i < m; ++i)
    v(i) = values[i];
  apf::DynamicVector b;
  apf::multiply(PT,v,b);
  apf::DynamicVector a;
  solveSVD(A,a,b);
  coeffs.setSize(n);
  for (int i=0; i < n; ++i)
    coeffs[i] = a[i];
}

double evalPolynomial(int order, apf::Vector3& point, apf::DynamicVector& coeffs)
{
  apf::DynamicVector terms;
  evalPolynomialTerms(order, point, terms);
  return coeffs * terms;
}

/** @brief run patch recovery
  */
void runSpr(Patch* p)
{
  Recovery* r = p->recovery;
  Samples samples;
  getSamples(p, &samples);
  int num_components = apf::countComponents(r->f_star);
  int num_nodes = r->mesh->getShape()->countNodesOn(r->mesh->getType(p->entity));
  apf::NewArray<double> values(samples.num_points);
  apf::NewArray<apf::Vector3> nodal_points(num_nodes);
  apf::NewArray<apf::NewArray<double> > recovered_values(num_nodes);
  for (int n = 0; n < num_nodes; ++n) {
    recovered_values[n].allocate(num_components);
    p->mesh->getPoint(p->entity, n, nodal_points[n]);
  }
  for (int i = 0; i < num_components; ++i) {
    for (int p = 0; p < samples.num_points; ++p)
      values[p] = samples.values[p][i];
    apf::DynamicVector coeffs;
    fitPolynomial(r->order, samples.num_points,
        samples.points, values, coeffs);
    for (int n = 0; n < num_nodes; ++n)
      recovered_values[n][i] = evalPolynomial(r->order, nodal_points[n], coeffs);
  }
  for (int n = 0; n < num_nodes; ++n)
    setComponents(r->f_star, p->entity, n, &(recovered_values[n][0]));
}

class PatchOp : public apf::CavityOp
{
public:
  PatchOp(Recovery* r):
    apf::CavityOp(r->mesh)
  {
    setupPatch(&patch, r);
  }
  virtual Outcome setEntity(apf::MeshEntity* e)
  {
    if (hasEntity(patch.recovery->f_star, e))
      return SKIP;
    startPatch(&patch, e);
    if ( ! buildPatch(&patch, this))
      return REQUEST;
    return OK;
  }
  virtual void apply()
  {
    runSpr(&patch);
  }
  Patch patch;
};

apf::Field* recoverField(apf::Field* f)
{
  Recovery recovery;
  setupRecovery(&recovery, f);
  PatchOp op(&recovery);
  for (int d = 0; d <= 3; ++d)
    if (recovery.mesh->getShape()->hasNodesIn(d))
      op.applyToDimension(d);
  return recovery.f_star;
}

}

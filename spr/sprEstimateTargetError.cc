/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

/* see Boussetta, Ramzy et al.
   "Adaptive remeshing based on a posteriori error estimation
   for forging simulation." Computer methods in applied mechanics
   and engineering 195.48 (2006): 6626-6645. */

#include "spr.h"

#include <PCU.h>
#include <apfMesh.h>
#include <apfShape.h>
#include <apfCavityOp.h>

#include <limits>
#include <cassert>

namespace spr {
namespace target {

struct Estimation {
  apf::Mesh* mesh;
  int integration_order;
  int recovered_order;
  apf::Field* eps;
  apf::Field* eps_star;
  size_t target_number;
  double size_factor;
  double alpha;
  double beta;
  apf::Field* elem_size;
  apf::Field* vtx_size;
};

static void setupEstimation(
    Estimation* e,
    apf::Field* eps,
    size_t target,
    double alpha,
    double beta)
{
  e->mesh = apf::getMesh(eps);
  e->integration_order = apf::getShape(eps)->getOrder();
  e->recovered_order = e->mesh->getShape()->getOrder();
  e->eps = eps;
  e->eps_star = 0;
  e->target_number = target;
  e->size_factor = 0;
  e->alpha = alpha;
  e->beta = beta;
  e->elem_size = 0;
  e->vtx_size = 0;
}

class ScalarIntegrator : public apf::Integrator
{
  public:
    ScalarIntegrator(int order):
      apf::Integrator(order),
      result(0)
    {
    }
    void parallelReduce()
    {
      PCU_Add_Doubles(&result,1);
    }
    double result;
};

class ElementError : public ScalarIntegrator
{
  public:
    ElementError(Estimation* e):
      ScalarIntegrator(e->integration_order),
      estimation(e),
      element(0),
      entity(0),
      sum(0),
      ip(0)
    {
      v1.setSize(apf::countComponents(e->eps));
      v2.setSize(apf::countComponents(e->eps_star));
    }
    void inElement(apf::MeshElement* meshElement)
    {
      element = apf::createElement(estimation->eps_star, meshElement);
      entity = apf::getMeshEntity(meshElement);
      sum = 0;
      ip = 0;
    }
    void outElement()
    {
      apf::destroyElement(element);
    }
    void atPoint(apf::Vector3 const& xi, double w, double dv)
    {
      apf::getComponents(estimation->eps, entity, ip, &v1[0]);
      apf::getComponents(element, xi, &v2[0]);
      apf::DynamicVector& diff = v1;
      diff -= v2;
      sum += (diff*diff)*w*dv;
      ++ip;
    }
    Estimation* estimation;
    apf::Element* element;
    apf::MeshEntity* entity;
    double sum;
    int ip;
    apf::DynamicVector v1,v2;
};

class GlobalErrorTerm : public ElementError
{
  public:
    GlobalErrorTerm(Estimation* e) : ElementError(e) {}
    void outElement()
    {
      ElementError::outElement();
      double d = estimation->mesh->getDimension();
      double p = estimation->recovered_order;
      result += pow(sqrt(sum), ((2*d)/(2*p+d)));
    }
};

static void computeSizeFactor(Estimation* e)
{
  GlobalErrorTerm globalError(e);
  globalError.process(e->mesh);
  int d = e->mesh->getDimension();
  double G = globalError.result;
  double N = e->target_number;
  e->size_factor = pow((G/N), (1.0/d));
}

static double getCurrentSize(apf::Mesh* m, apf::MeshEntity* e)
{
  double h=0;
  apf::Downward edges;
  int ne = m->getDownward(e,1,edges);
  for (int i=0; i < ne; ++i)
    h = std::max(h, apf::measure(m, edges[i]));
  return h;
}

static double getNewSize(Estimation* e, apf::MeshEntity* ent)
{
  ElementError elemError(e);
  apf::MeshElement* elem = apf::createMeshElement(e->mesh, ent);
  elemError.process(elem);
  int p = e->recovered_order;
  int d = e->mesh->getDimension();
  double h = getCurrentSize(e->mesh, ent);
  double theta_e = sqrt(elemError.sum);
  double r = pow(theta_e, ((-2.0)/(2.0*p+d)));
  double h_new = e->size_factor * r * h;
  if (h_new < e->alpha*h) h_new = e->alpha*h;
  if (h_new > e->beta*h) h_new = e->beta*h;
  return h_new;
}

static void getElemSizeField(Estimation* e)
{
  apf::Field* esize = apf::createStepField(e->mesh,"esize",apf::SCALAR);
  int d = e->mesh->getDimension();
  apf::MeshEntity* elem;
  apf::MeshIterator* elems = e->mesh->begin(d);
  while ((elem = e->mesh->iterate(elems))) {
    double h = getNewSize(e, elem);
    apf::setScalar(esize, elem, 0, h);
  }
  e->mesh->end(elems);
  e->elem_size = esize;
}

static void avgToVtx(apf::Field* ef, apf::Field* vf, apf::MeshEntity* ent)
{
  apf::Mesh* m = apf::getMesh(ef);
  apf::Adjacent elems;
  m->getAdjacent(ent, m->getDimension(), elems);
  double s = 0;
  for (size_t i=0; i < elems.getSize(); ++i)
    s += apf::getScalar(ef, elems[i], 0);
  s /= elems.getSize();
  apf::setScalar(vf, ent, 0, s);
}

class AverageOp : public apf::CavityOp
{
  public:
    AverageOp(Estimation* est):
      apf::CavityOp(est->mesh),
      estimation(est),
      entity(0)
    {
    }
    virtual Outcome setEntity(apf::MeshEntity* e)
    {
      entity = e;
      if (apf::hasEntity(estimation->vtx_size, entity))
        return SKIP;
      if ( ! requestLocality(&entity, 1))
        return REQUEST;
      return OK;
    }
    virtual void apply()
    {
      avgToVtx(
          estimation->elem_size,
          estimation->vtx_size,
          entity);
    }
    Estimation* estimation;
    apf::MeshEntity* entity;
};

static void averageSizeField(Estimation* e)
{
  e->vtx_size =
    apf::createLagrangeField(e->mesh, "size", apf::SCALAR, 1);
  AverageOp op(e);
  op.applyToDimension(0);
}

static void estimateError(Estimation* e)
{
  e->eps_star = recoverField(e->eps);
  computeSizeFactor(e);
  getElemSizeField(e);
  apf::destroyField(e->eps_star);
  averageSizeField(e);
  apf::destroyField(e->elem_size);
}

}

apf::Field* getTargetSPRSizeField(
    apf::Field* eps,
    size_t target,
    double alpha,
    double beta)
{
  double t0 = PCU_Time();
  assert(target > 0);
  assert(alpha < beta);
  target::Estimation e;
  target::setupEstimation(&e, eps, target, alpha, beta);
  target::estimateError(&e);
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
    fprintf(stderr, "SPR (target): error estimated in %f seconds\n",t1-t0);
  return e.vtx_size;
}

}

/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <PCU.h>

#include "spr.h"

#include <apfMesh.h>
#include <apfShape.h>
#include <apfCavityOp.h>

#include <limits>

namespace spr {

/* common base for Scalar Integrator. */
class SInt : public apf::Integrator
{
  public:
    SInt(int order):
      apf::Integrator(order),r(0)
    {}
    void parallelReduce()
    {
      PCU_Add_Doubles(&r,1);
    }
    void reset() {r=0;}
    double r;
};

struct Estimation {
  apf::Mesh* mesh;
  /* the maximum polynomial order that can be
     integrated exactly using the input integration points.
     not necessarily equal to recovered_order, sometimes
     users give more points than strictly necessary */
  int integration_order;
  /* the polynomial order of the recovered field */
  int recovered_order;
  /* the input field consisting of values of a key
     quantity at integration points (stress or strain, for example).
     this field is related to integration_order */
  apf::Field* eps;
  /* the recovered field, consisting of nodal values
     of the quantity of interest.
     this uses a Lagrange basis of recovered_order */
  apf::Field* eps_star;
  /* the acceptable margin of error, expressed as a factor
     greater than zero.
     setting this equal to zero would request zero element
     size everywhere, i.e. infinite refinement.
     increasing it scales up the desired element size
     throughout.
     basically, this is a (nonlinear) scaling factor on the resulting
     size field */
  double tolerance;
  /* the uniform linear scaling factor derived from the tolerance
     and integrals of the recovered field over the mesh.
     desired element size = current size * current error * size_factor */
  double size_factor;
  /* a temporary field storing desired sizes at elements */
  apf::Field* element_size;
  /* the resulting size field, recovered from the element_size field
     (using a local average recovery method much weaker than SPR) */
  apf::Field* size;
};

/* useful for initializing values to quickly
   detect "uninitialized" value bugs */
static double getNaN()
{
  return std::numeric_limits<double>::quiet_NaN();
}

static void setupEstimation(Estimation* e, apf::Field* eps, double tolerance)
{
  /* note that getOrder being used to convey this meaning
     is a bit of a hack, but looks decent if you don't
     try to define the FieldShape API too rigorously */
  e->integration_order = apf::getShape(eps)->getOrder();
  e->mesh = apf::getMesh(eps);
  /* so far recovery order is directly tied to the
     mesh's coordinate field order, coordinate this
     with field recovery code */
  e->recovered_order = e->mesh->getShape()->getOrder();
  e->eps = eps;
  e->tolerance = tolerance;
  e->size_factor = getNaN();
  e->element_size = 0;
  e->size = 0;
}

/* computes $\|f\|^2$ */
class SelfProduct : public SInt
{
  public:
    SelfProduct(Estimation* e):
      SInt(e->integration_order), estimation(e)
    {
      v.setSize(apf::countComponents(e->eps_star));
    }
    void inElement(apf::MeshElement* meshElement)
    {
      element = apf::createElement(estimation->eps_star, meshElement);
    }
    void outElement()
    {
      apf::destroyElement(element);
    }
    void atPoint(apf::Vector3 const& p, double w, double dV)
    {
      apf::getComponents(element, p, &v[0]);
      r += (v * v) * w * dV;
    }
  private:
    Estimation* estimation;
    apf::Element* element;
    apf::DynamicVector v;
};

/* computes the integral over the element of the
   sum of the squared differences between the
   original and recovered fields */
class ElementError : public SInt
{
  public:
    ElementError(Estimation* e):
      SInt(e->integration_order)
    {
      estimation = e;
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
    void atPoint(apf::Vector3 const& xi, double w, double dV)
    {
      apf::getComponents(estimation->eps, entity, ip, &v1[0]);
      apf::getComponents(element, xi, &v2[0]);
      apf::DynamicVector& diff = v1;
      diff -= v2;
      sum += (diff * diff) * w * dV;
      ++ip;
    }
    Estimation* estimation;
    apf::Element* element;
    apf::MeshEntity* entity;
    double sum; //element result
    int ip; //integration point counter
    apf::DynamicVector v1, v2;
};

/* computes the $\sum_{i=1}^n \|e_\epsilon\|^{\frac{2d}{2p+d}}$ term. */
class Error : public ElementError
{
  public:
    Error(Estimation* e):
      ElementError(e)
    {
    }
    void outElement()
    {
      ElementError::outElement();
      double d = estimation->mesh->getDimension();
      double p = estimation->recovered_order;
      r += pow(sqrt(sum), ((2 * d) / (2 * p + d)));
    }
};

/* computes the $\|e_\epsilon\|^{-\frac{2}{2p+d}}_e$ term
   (over one element only) */
class Error2 : public ElementError
{
  public:
    Error2(Estimation* e):
      ElementError(e)
    {
    }
    void outElement()
    {
      ElementError::outElement();
      double p = estimation->recovered_order;
      double d = estimation->mesh->getDimension();
      r = pow(sqrt(sum), -(2 / (2 * p + d)));
    }
};

static void computeSizeFactor(Estimation* e)
{
  SelfProduct epsStarNormIntegrator(e);
  epsStarNormIntegrator.process(e->mesh);
  double epsStarNorm = sqrt(epsStarNormIntegrator.r);
  Error errorIntegrator(e);
  errorIntegrator.process(e->mesh);
  double a = e->tolerance * e->tolerance *
             epsStarNorm * epsStarNorm;
  double b = a / errorIntegrator.r;
  double p = e->recovered_order;
  e->size_factor = pow(b, 1.0 / (2.0 * p));
}

static double getEdgeLength(apf::Mesh* m, apf::MeshEntity* e)
{
  apf::MeshElement* element = apf::createMeshElement(m,e);
  double h = measure(element);
  apf::destroyMeshElement(element);
  return h;
}

static double getCurrentSize(apf::Mesh* m, apf::MeshEntity* e)
{
  /* right now maximum edge length is the formula... */
  double h = 0;
  apf::Downward edges;
  int ne = m->getDownward(e,1,edges);
  for (int i=0; i < ne; ++i)
    h = std::max(h, getEdgeLength(m, edges[i]));
  return h;
}

static double getDesiredSize(Estimation* e, apf::MeshEntity* entity)
{
  Error2 errorNormIntegrator(e);
  apf::MeshElement* element = apf::createMeshElement(e->mesh, entity);
  errorNormIntegrator.process(element);
  double errorNorm = errorNormIntegrator.r;
  apf::destroyMeshElement(element);
  double h = getCurrentSize(e->mesh, entity);
  return h * errorNorm * e->size_factor;
}

static void getElementSizeField(Estimation* e)
{
  apf::Field* eSize = apf::createStepField(e->mesh, "esize", apf::SCALAR);
  int d = e->mesh->getDimension();
  apf::MeshEntity* entity;
  apf::MeshIterator* elements = e->mesh->begin(d);
  while ((entity = e->mesh->iterate(elements))) {
    double h = getDesiredSize(e, entity);
    apf::setScalar(eSize, entity, 0, h);
  }
  e->mesh->end(elements);
  e->element_size = eSize;
}

/* note that this only works when there is
   one node per entity at most. */
void averageToEntity(apf::Field* ef, apf::Field* vf, apf::MeshEntity* ent)
{
  apf::Mesh* m = apf::getMesh(ef);
  apf::Adjacent elements;
  m->getAdjacent(ent, m->getDimension(), elements);
  double s=0;
  for (std::size_t i=0; i < elements.getSize(); ++i)
    s += apf::getScalar(ef, elements[i], 0);
  s /= elements.getSize();
  apf::setScalar(vf, ent, 0, s);
}

class AverageOp : public apf::CavityOp
{
  public:
    AverageOp(Estimation* e):
      apf::CavityOp(e->mesh)
    {
      estimation = e;
    }
    virtual Outcome setEntity(apf::MeshEntity* e)
    {
      entity = e;
      if (apf::hasEntity(estimation->size, entity))
        return SKIP;
      if ( ! requestLocality(&entity, 1))
        return REQUEST;
      return OK;
    }
    virtual void apply()
    {
      averageToEntity(estimation->element_size,
          estimation->size, entity);
    }
    Estimation* estimation;
    apf::MeshEntity* entity;
};

void averageSizeField(Estimation* e)
{
  e->size = apf::createFieldOn(e->mesh, "size", apf::SCALAR);
  AverageOp op(e);
  for (int d = 0; d <= e->mesh->getDimension(); ++d)
    if (e->mesh->getShape()->hasNodesIn(d))
      op.applyToDimension(d);
}

static void estimateError(Estimation* e)
{
  e->eps_star = recoverField(e->eps);
  computeSizeFactor(e);
  getElementSizeField(e);
  apf::destroyField(e->eps_star);
  averageSizeField(e);
  apf::destroyField(e->element_size);
}

apf::Field* getSPRSizeField(apf::Field* eps, double adaptRatio)
{
  double t0 = PCU_Time();
  Estimation e;
  setupEstimation(&e, eps, adaptRatio);
  estimateError(&e);
  double t1 = PCU_Time();
  if (!PCU_Comm_Self())
    fprintf(stderr,"SPR: error estimated in %f seconds\n",t1-t0);
  return e.size;
}

}

/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "spr.h"
#include "apfMesh.h"
#include "apfShape.h"
#include "apfCavityOp.h"
#include <PCU.h>
#include <algorithm>

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

/* computes \|f\|^2 */
class SelfProduct : public SInt
{
  public:
    SelfProduct(apf::Field* f_i,int order):
      SInt(order),f(f_i),e(0) {}
    void inElement(apf::MeshElement* meshElement)
    {
      e = apf::createElement(f,meshElement);
    }
    void outElement()
    {
      apf::destroyElement(e);
    }
    void atPoint(apf::Vector3 const& p, double w, double dV)
    {
      apf::Matrix3x3 v;
      apf::getMatrix(e,p,v);
      r += apf::getInnerProduct(v,v)*w*dV;
    }
  private:
    apf::Field* f;
    apf::Element* e;
};

/* computes the $\sum_{i=1}^n \|e_\epsilon\|^{\frac{2d}{2p+d}}$ term. */
class Error : public SInt
{
  public:
    Error(apf::Field* e, apf::Field* es, int order):
      SInt(order),eps(e),eps_star(es),e(0),p(order)
    {
      apf::Mesh* m = apf::getMesh(eps);
      d = m->getDimension();
    }
    void inElement(apf::MeshElement* meshElement)
    {
      e = apf::createElement(eps_star,meshElement);
      er = 0;
      ip = 0;
    }
    void outElement()
    {
      r += pow(sqrt(er),((2*d)/(2*p + d)));
      apf::destroyElement(e);
    }
    void atPoint(apf::Vector3 const& xi, double w, double dV)
    {
      apf::Matrix3x3 v1,v2,diff;
      apf::getMatrix(eps,apf::getMeshEntity(apf::getMeshElement(e)),ip,v1);
      apf::getMatrix(e,xi,v2);
      diff = v1-v2;
      er += apf::getInnerProduct(diff,diff)*w*dV;
      ++ip;
    }
  private:
    apf::Field* eps; //IP sampled gradient
    apf::Field* eps_star; //recovered gradient
    apf::Element* e; //recovered gradient element
    double er; //element result
    double p,d; //polynomial order and element dimension
    int ip; //integration point counter
};

/* computes the $\|e_\epsilon\|^{-\frac{2}{2p+d}}_e$ term
   (over elements only) */
class Error2 : public SInt
{
  public:
    Error2(apf::Field* e, apf::Field* es, int order):
      SInt(order),eps(e),eps_star(es),e(0),p(order)
    {
      apf::Mesh* m = apf::getMesh(eps);
      d = m->getDimension();
    }
    void inElement(apf::MeshElement* meshElement)
    {
      e = apf::createElement(eps_star,meshElement);
      ip = 0;
    }
    void outElement()
    {
      r = pow(sqrt(r),-(2/(2*p + d)));
      apf::destroyElement(e);
    }
    void atPoint(apf::Vector3 const& p, double w, double dV)
    {
      apf::Matrix3x3 v1,v2,diff;
      apf::getMatrix(eps,apf::getMeshEntity(apf::getMeshElement(e)),ip,v1);
      apf::getMatrix(e,p,v2);
      diff = v1-v2;
      r += apf::getInnerProduct(diff,diff)*w*dV;
      ++ip;
    }
  private:
    apf::Field* eps; //IP sampled gradient
    apf::Field* eps_star; //recovered gradient
    apf::Element* e; //recovered gradient element
    double p,d; //polynomial order and element dimension
    int ip; //integration point counter
};

double computeSizeFactor(apf::Field* eps,
                         apf::Field* eps_star,
                         double adaptRatio)
{
  apf::Mesh* mesh = apf::getMesh(eps);
  double p = mesh->getShape()->getOrder();
  SelfProduct epsStarNormIntegrator(eps_star,p);
  epsStarNormIntegrator.process(mesh);
  double epsStarNorm = sqrt(epsStarNormIntegrator.r);
  Error errorIntegrator(eps,eps_star,p);
  errorIntegrator.process(mesh);
  return pow((((adaptRatio*adaptRatio)*
               (epsStarNorm*epsStarNorm))/
              errorIntegrator.r),
             1.0/(2.0*p));
}

double getEdgeLength(apf::Mesh* m, apf::MeshEntity* e)
{
  apf::MeshElement* element = apf::createMeshElement(m,e);
  double h = measure(element);
  apf::destroyMeshElement(element);
  return h;
}

double getCurrentSize(apf::Mesh* m, apf::MeshEntity* e)
{
  /* right now maximum edge length is the formula... */
  double h = 0;
  apf::Downward edges;
  int ne = m->getDownward(e,1,edges);
  for (int i=0; i < ne; ++i)
    h = std::max(h,getEdgeLength(m,edges[i]));
  return h;
}

double getDesiredSize(apf::Mesh* m,
                      apf::MeshEntity* e,
                      apf::Field* eps, apf::Field* eps_star,
                      double sizeFactor)
{
  int p = m->getShape()->getOrder();
  Error2 errorNormIntegrator(eps,eps_star,p);
  apf::MeshElement* element = apf::createMeshElement(m,e);
  errorNormIntegrator.process(element);
  double errorNorm = errorNormIntegrator.r;
  apf::destroyMeshElement(element);
  double h = getCurrentSize(m,e);
  return h*errorNorm*sizeFactor;
}

apf::Field* getElementSizeField(apf::Field* eps,
                           apf::Field* eps_star,
                           double sizeFactor)
{
  apf::Mesh* mesh = apf::getMesh(eps);
  apf::Field* eSize = apf::createStepField(mesh,"esize",apf::SCALAR);
  int d = mesh->getDimension();
  apf::MeshEntity* entity;
  apf::MeshIterator* elements = mesh->begin(d);
  while ((entity = mesh->iterate(elements)))
  {
    double h = getDesiredSize(mesh,entity,eps,eps_star,sizeFactor);
    apf::setScalar(eSize,entity,0,h);
  }
  mesh->end(elements);
  return eSize;
}

void averageToVertices(apf::Field* ef, apf::Field* vf, apf::MeshEntity* v)
{
  apf::Mesh* m = apf::getMesh(ef);
  apf::Adjacent elements;
  m->getAdjacent(v,m->getDimension(),elements);
  double s=0;
  for (std::size_t i=0; i < elements.getSize(); ++i)
    s += apf::getScalar(ef,elements[i],0);
  s /= elements.getSize();
  apf::setScalar(vf,v,0,s);
}

class ElementsToVertex : public apf::CavityOp
{
  public:
    ElementsToVertex(apf::Field* e, apf::Field* v):
      apf::CavityOp(apf::getMesh(e)),
      elementField(e),
      vertexField(v)
    {}
    virtual Outcome setEntity(apf::MeshEntity* e)
    {
      vertex = e;
      if (hasEntity(vertexField,vertex))
        return SKIP;
      if ( ! requestLocality(&vertex,1))
        return REQUEST;
      return OK;
    }
    virtual void apply()
    {
      averageToVertices(elementField,vertexField,vertex);
    }
    apf::Field* elementField;
    apf::Field* vertexField;
    apf::MeshEntity* vertex;
};

apf::Field* elementToVertexField(apf::Field* eSize)
{
  apf::Mesh* m = apf::getMesh(eSize);
  apf::Field* sizeField = apf::createFieldOn(m,"size",apf::SCALAR);
  ElementsToVertex op(eSize,sizeField);
  op.applyToDimension(0);
  /* averaging vtx size field to edge node for now
     ElementsToVertex should really probably be
     extended to ElementsToNode */
  if (m->getShape()->countNodesOn(1) != 0)
  {
    apf::MeshIterator* edges = m->begin(1);
    apf::MeshEntity* edge;
    while ((edge = m->iterate(edges)))
    {
      apf::MeshEntity* v[2];
      m->getDownward(edge,0,v);
      double s[2];
      s[0] = apf::getScalar(sizeField,v[0],0);
      s[1] = apf::getScalar(sizeField,v[1],0);
      double nodeSize = 0.5*(s[0]+s[1]);
      apf::setScalar(sizeField,edge,0,nodeSize);
    }
    m->end(edges);
  }
  return sizeField;
}

apf::Field* getSPRSizeField(apf::Field* eps, double adaptRatio)
{
  apf::Field* eps_star = recoverField(eps);
  double sizeFactor = computeSizeFactor(eps,eps_star,adaptRatio);
  apf::Field* eSize = getElementSizeField(eps,eps_star,sizeFactor);
  apf::destroyField(eps_star);
  apf::Field* sizeField;
  sizeField = elementToVertexField(eSize);
  apf::destroyField(eSize);
  return sizeField;
}

}

/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <algorithm>
#include <PCU.h>
#include "spr.h"
#include "apfMesh.h"
#include "apfShape.h"
#include "apfCavityOp.h"

namespace apf {

/* common base for Scalar Integrator. */
class SInt : public Integrator
{
  public:
    SInt(int order):
      Integrator(order),r(0)
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
    SelfProduct(Field* f_i,int order):
      SInt(order),f(f_i),e(0) {}
    void inElement(MeshElement* meshElement)
    {
      e = createElement(f,meshElement);
    }
    void outElement()
    {
      destroyElement(e);
    }
    void atPoint(Vector3 const& p, double w, double dV)
    {
      Matrix3x3 v;
      getMatrix(e,p,v);
      r += getInnerProduct(v,v)*w*dV;
    }
  private:
    Field* f;
    Element* e;
};

/* computes the $\sum_{i=1}^n \|e_\epsilon\|^{\frac{2d}{2p+d}}$ term. */
class Error : public SInt
{
  public:
    Error(Field* e, Field* es, int order):
      SInt(order),eps(e),eps_star(es),e(0),p(order)
    {
      Mesh* m = getMesh(eps);
      d = m->getDimension();
    }
    void inElement(MeshElement* meshElement)
    {
      e = createElement(eps_star,meshElement);
      er = 0;
      ip = 0;
    }
    void outElement()
    {
      r += pow(sqrt(er),((2*d)/(2*p + d)));
      destroyElement(e);
    }
    void atPoint(Vector3 const& xi, double w, double dV)
    {
      Matrix3x3 v1,v2,diff;
      getMatrix(eps,getMeshEntity(getMeshElement(e)),ip,v1);
      getMatrix(e,xi,v2);
      diff = v1-v2;
      er += getInnerProduct(diff,diff)*w*dV;
      ++ip;
    }
  private:
    Field* eps; //IP sampled gradient
    Field* eps_star; //recovered gradient
    Element* e; //recovered gradient element
    double er; //element result
    double p,d; //polynomial order and element dimension
    int ip; //integration point counter
};

/* computes the $\|e_\epsilon\|^{-\frac{2}{2p+d}}_e$ term
   (over elements only) */
class Error2 : public SInt
{
  public:
    Error2(Field* e, Field* es, int order):
      SInt(order),eps(e),eps_star(es),e(0),p(order)
    {
      Mesh* m = getMesh(eps);
      d = m->getDimension();
    }
    void inElement(MeshElement* meshElement)
    {
      e = createElement(eps_star,meshElement);
      ip = 0;
    }
    void outElement()
    {
      r = pow(sqrt(r),-(2/(2*p + d)));
      destroyElement(e);
    }
    void atPoint(Vector3 const& p, double w, double dV)
    {
      Matrix3x3 v1,v2,diff;
      getMatrix(eps,getMeshEntity(getMeshElement(e)),ip,v1);
      getMatrix(e,p,v2);
      diff = v1-v2;
      r += getInnerProduct(diff,diff)*w*dV;
      ++ip;
    }
  private:
    Field* eps; //IP sampled gradient
    Field* eps_star; //recovered gradient
    Element* e; //recovered gradient element
    double p,d; //polynomial order and element dimension
    int ip; //integration point counter
};

double computeSizeFactor(Field* eps,
                         Field* eps_star,
                         double adaptRatio)
{
  Mesh* mesh = getMesh(eps);
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

double getEdgeLength(Mesh* m, MeshEntity* e)
{
  MeshElement* element = createMeshElement(m,e);
  double h = measure(element);
  destroyMeshElement(element);
  return h;
}

double getCurrentSize(Mesh* m, MeshEntity* e)
{
  /* right now maximum edge length is the formula... */
  double h = 0;
  Downward edges;
  int ne = m->getDownward(e,1,edges);
  for (int i=0; i < ne; ++i)
    h = std::max(h,getEdgeLength(m,edges[i]));
  return h;
}

double getDesiredSize(Mesh* m,
                      MeshEntity* e,
                      Field* eps, Field* eps_star,
                      double sizeFactor)
{
  int p = m->getShape()->getOrder();
  Error2 errorNormIntegrator(eps,eps_star,p);
  MeshElement* element = createMeshElement(m,e);
  errorNormIntegrator.process(element);
  double errorNorm = errorNormIntegrator.r;
  destroyMeshElement(element);
  double h = getCurrentSize(m,e);
  return h*errorNorm*sizeFactor;
}

Field* getElementSizeField(Field* eps,
                           Field* eps_star,
                           double sizeFactor)
{
  Mesh* mesh = getMesh(eps);
  Field* eSize = createStepField(mesh,"esize",SCALAR);
  int d = mesh->getDimension();
  MeshEntity* entity;
  MeshIterator* elements = mesh->begin(d);
  while ((entity = mesh->iterate(elements)))
  {
    double h = getDesiredSize(mesh,entity,eps,eps_star,sizeFactor);
    setScalar(eSize,entity,0,h);
  }
  mesh->end(elements);
  return eSize;
}

void averageToVertices(Field* ef, Field* vf, MeshEntity* v)
{
  Mesh* m = getMesh(ef);
  Adjacent elements;
  m->getAdjacent(v,m->getDimension(),elements);
  double s=0;
  for (std::size_t i=0; i < elements.getSize(); ++i)
    s += getScalar(ef,elements[i],0);
  s /= elements.getSize();
  setScalar(vf,v,0,s);
}

class ElementsToVertex : public CavityOp
{
  public:
    ElementsToVertex(Field* e, Field* v):
      CavityOp(getMesh(e)),
      elementField(e),
      vertexField(v)
    {}
    virtual Outcome setEntity(MeshEntity* e)
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
    Field* elementField;
    Field* vertexField;
    MeshEntity* vertex;
};

Field* elementToVertexField(Field* eSize)
{
  Mesh* m = getMesh(eSize);
  Field* sizeField = createFieldOn(m,"size",SCALAR);
  ElementsToVertex op(eSize,sizeField);
  op.applyToDimension(0);
  /* averaging vtx size field to edge node for now
     ElementsToVertex should really probably be
     extended to ElementsToNode */
  if (m->getShape()->countNodesOn(1) != 0)
  {
    MeshIterator* edges = m->begin(1);
    MeshEntity* edge;
    while ((edge = m->iterate(edges)))
    {
      MeshEntity* v[2];
      m->getDownward(edge,0,v);
      double s[2];
      s[0] = getScalar(sizeField,v[0],0);
      s[1] = getScalar(sizeField,v[1],0);
      double nodeSize = 0.5*(s[0]+s[1]);
      setScalar(sizeField,edge,0,nodeSize);
    }
    m->end(edges);
  }
  return sizeField;
}

Field* getSPRSizeField(Field* eps, double adaptRatio)
{
  Field* eps_star = recoverField(eps);
  double sizeFactor = computeSizeFactor(eps,eps_star,adaptRatio);
  Field* eSize = getElementSizeField(eps,eps_star,sizeFactor);
  destroyField(eps_star);
  Field* sizeField;
  sizeField = elementToVertexField(eSize);
  destroyField(eSize);
  return sizeField;
}

}//namespace apf

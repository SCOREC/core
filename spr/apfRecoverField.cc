/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apf.h"
#include "apfMesh.h"
#include "apfShape.h"
#include "apfSVDecomp.h"
#include "apfCavityOp.h"
#include <set>

namespace apf {

Field* getGradIPField(Field* f, const char* name, int order)
{
  Mesh* m = getMesh(f);
  int vt = getValueType(f);
  assert(vt == SCALAR || vt == VECTOR);
  Field* ipField = createIPField(m,name,vt+1,order);
  MeshIterator* it = m->begin(m->getDimension());
  MeshEntity* e;
  while ((e = m->iterate(it)))
  {
    MeshElement* me = createMeshElement(m,e);
    Element* fe = createElement(f,me);
    int np = countIntPoints(me,order);
    for (int p=0; p < np; ++p)
    {
      Vector3 xi;
      getIntPoint(me,order,p,xi);
      if (vt == SCALAR)
      {
        Vector3 value;
        getGrad(fe,xi,value);
        setVector(ipField,e,p,value);
      }
      else
      {
        Matrix3x3 value;
        getVectorGrad(fe,xi,value);
        setMatrix(ipField,e,p,value);
      }
    }
    destroyElement(fe);
    destroyMeshElement(me);
  }
  m->end(it);
  return ipField;
}

struct SprPoints {
  SprPoints():count(0) {}
  void setCount(int n)
  {
    count = n;
    points.allocate(n);
    values.allocate(n);
  }
  int count;
  NewArray<Vector3> points;
  NewArray<Matrix3x3> values;
};

void evalPolynomialTerms(int order,
                         Vector3 const& point,
                         DynamicVector& polynomial)
{
  Vector3 const& x = point;
  if (order == 1)
  {
    polynomial.setSize(4);
    polynomial(0) = 1.0;
    polynomial(1) = x[0];
    polynomial(2) = x[1];
    polynomial(3) = x[2];
  }
  else if (order == 2)
  {
    polynomial.setSize(10);
    polynomial(0) = 1.0;
    polynomial(1) = x[0];
    polynomial(2) = x[1];
    polynomial(3) = x[2];
    polynomial(4) = x[0]*x[1];
    polynomial(5) = x[1]*x[2];
    polynomial(6) = x[2]*x[0];
    polynomial(7) = x[0]*x[0];
    polynomial(8) = x[1]*x[1];
    polynomial(9) = x[2]*x[2];
  }
  else
    apf::fail("SPR: invalid polynomial order");
}

void fitPolynomial(int count,
                   int order,
                   NewArray<Vector3> const& points,
                   NewArray<double> const& values,
                   DynamicVector& coefficients)
{
  DynamicVector v(count);
  for (int i=0; i < count; ++i)
    v(i) = values[i];
  int m = (order+1)*(order+2)*(order+3)/6;
  DynamicMatrix P(count,m);
  DynamicVector p;
  for (int i = 0; i < count; ++i)
  {
    evalPolynomialTerms(order,points[i],p);
    P.setRow(i,p);
  }
  DynamicMatrix PT;
  transpose(P,PT);
  DynamicMatrix A;
  multiply(PT,P,A);
  DynamicVector b;
  multiply(PT,v,b);
  DynamicVector a;
  /* HELP: both these values are random and we need
     to get good values with supporting evidence */
  solveSVD(A,a,b);
  coefficients.setSize(m);
  for (int i=0; i < m; ++i)
    coefficients[i] = a[i];
}

void getSprPoints(Field* eps,
                  std::set<MeshEntity*>& elements,
                  SprPoints& sprPoints)
{
  Mesh* mesh = getMesh(eps);
  std::set<MeshEntity*>::iterator elem = elements.begin();
  MeshElement* me = createMeshElement(mesh,*elem);
  int order = getOrder(me);
  int pointCount = countIntPoints(me,order);
  destroyMeshElement(me);
  /* assumes const #IP across elements */
  sprPoints.setCount(elements.size()*pointCount);
  std::size_t p = 0;
  for (std::set<MeshEntity*>::iterator it = elements.begin();
      it != elements.end(); ++it)
  {
    MeshElement* meshElement = createMeshElement(mesh,*it);
    for (int l=0; l < pointCount; ++l)
    {
      Vector3 param;
      getIntPoint(meshElement,order,l,param);
      mapLocalToGlobal(meshElement,param,sprPoints.points[p]);
      getMatrix(eps,*it,l,sprPoints.values[p]);
      ++p;
    }
    destroyMeshElement(meshElement);
  }
}

void runSpr(Field* eps,
            Field* eps_star,
            std::set<MeshEntity*>& elements,
            MeshEntity* entity)
{
  SprPoints sprPoints;
  getSprPoints(eps,elements,sprPoints);
  NewArray<double> componentValues(sprPoints.count);
  Vector3 point;
  getMesh(eps)->getPoint(entity,0,point);
  Matrix3x3 value;
  Mesh* mesh = getMesh(eps);
  int order = mesh->getShape()->getOrder();
  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
    {
      for (int p=0; p < sprPoints.count; ++p)
        componentValues[p] = sprPoints.values[p][i][j];
      DynamicVector polynomial;
      DynamicVector coefficients;
      fitPolynomial(sprPoints.count,
                    order,
                    sprPoints.points,
                    componentValues,
                    coefficients);
      evalPolynomialTerms(order,point,polynomial);
      value[i][j] = coefficients*polynomial;
    }
  /* Need to set a matrix value to each node on a mesh
     entity in order to generalize. OK for 1st and 2nd
     order for now */
  setMatrix(eps_star,entity,0,value);
}

class PatchOp : public CavityOp
{
  public:
    PatchOp(Field* e, Field* es):
      CavityOp(getMesh(e)),
      eps(e),
      eps_star(es),
      pointCount(0)
    {
    }
    void reset(MeshEntity* e)
    {
      pointCount = 0;
      elements.clear();
      entity = e;
    }
    virtual Outcome setEntity(MeshEntity* e)
    {
      if (hasEntity(eps_star,e))
        return SKIP;
      reset(e);
      if ( ! buildPatch())
        return REQUEST;
      return OK;
    }
    virtual void apply()
    {
      runSpr(eps,eps_star,elements,entity);
    }
    bool buildPatch()
    {
      if (!getInitialPatch()) return false;
      if (!expandAsNecessary()) return false;
      return true;
    }
    bool getInitialPatch()
    {
      if ( ! requestLocality(&entity,1)) return false;
      DynamicArray<MeshEntity*> adjacent;
      mesh->getAdjacent(entity,mesh->getDimension(),adjacent);
      add(adjacent);
      return true;
    }
    bool expandAsNecessary()
    {
      const int desiredCount = 4;
      if (pointCount >= desiredCount)
        return true;
      std::set<MeshEntity*> oldSet = elements;
      int d = mesh->getDimension();
      for (int shared_dim=d-1; shared_dim >= 0; --shared_dim)
      {
        if (!addElementsThatShare(shared_dim,oldSet))
          return false;
        if (pointCount >= desiredCount)
          return true;
      }
      bool hope = elements.size() > oldSet.size();
      if (hope)
        return expandAsNecessary();
      else
      {
        fail("SPR Patch Builder: all hope is lost.");
        return false;
      }
    }
    bool addElementsThatShare(int dimension,
                              std::set<MeshEntity*>& oldSet)
    {
      std::set<MeshEntity*> bridges;
      for (std::set<MeshEntity*>::iterator it = oldSet.begin();
          it != oldSet.end(); ++it)
      {
        Downward down;
        int nd = mesh->getDownward(*it,dimension,down);
        for (int i=0; i < nd; ++i)
          bridges.insert(down[i]);
      }
      std::vector<MeshEntity*> bridgeArray(bridges.begin(),bridges.end());
      bridges.clear();
      if ( ! requestLocality(&(bridgeArray[0]),bridgeArray.size()))
        return false;
      for (size_t i=0; i < bridgeArray.size(); ++i)
      {
        Adjacent candidates;
        mesh->getAdjacent(bridgeArray[i],mesh->getDimension(),candidates);
        add(candidates);
      }
      return true;
    }
    void add(DynamicArray<MeshEntity*>& es)
    {
      for (std::size_t i=0; i < es.getSize(); ++i)
        add(es[i]);
    }
    void add(MeshEntity* e)
    {
      if (elements.count(e))
        return;
      elements.insert(e);
      MeshElement* element = createMeshElement(mesh,e);
      pointCount += countIntPoints(element,1);
      destroyMeshElement(element);
    }
    Field* eps;
    Field* eps_star;
    int pointCount;
    std::set<MeshEntity*> elements;
    MeshEntity* entity;
};

Field* recoverField(Field* eps)
{
  std::string name = "spr_";
  name += getName(eps);
  Mesh* m = getMesh(eps);
  Field* eps_star = createFieldOn(m,name.c_str(),MATRIX);
  PatchOp op(eps,eps_star);
  for (int i=0; i <= 3; ++i)
  {
    if (m->getShape()->countNodesOn(i) != 0)
      op.applyToDimension(i);
  }
  return eps_star;
}

}

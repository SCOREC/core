/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apf.h"
#include "apfMesh.h"
#include "apfField.h"
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

/** @brief least squares solution of polynomial coefficients
  * @param order (In) order of polynomial to be fit to data
  * @param num_points (In) number of integration points in patch
  * @param points (In) coordinates of integration points in patch
  * @param values (In) data value associated with integration points
  * @param coeffs (Out) coefficients of fit polynomial
  */
void evalPolynomialCoeffs(int order,
                          int num_points,
                          NewArray<Vector3> const& points,
                          NewArray<double> const& values,
                          DynamicVector& coeffs)
{
  int m = num_points;
  int n = (order+1)*(order+2)*(order+3)/6;
  DynamicMatrix P(m,n);
  DynamicVector p;
  for (int i = 0; i < m; ++i)
  {
    evalPolynomialTerms(order,points[i],p);
    P.setRow(i,p);
  }
  DynamicMatrix PT;
  transpose(P,PT);
  DynamicMatrix A;
  multiply(PT,P,A);
  DynamicVector v(m);
  for (int i=0; i < m; ++i)
    v(i) = values[i];
  DynamicVector b;
  multiply(PT,v,b);
  DynamicVector a;
  solveSVD(A,a,b);
  coeffs.setSize(n);
  for (int i=0; i < n; ++i)
    coeffs[i] = a[i];
}

/** @brief get spr point data from element patch
  * @param f (In) integration point field
  * @param elements (In) elements in patch
  * @param spr_points (Out) coord and data at points
  * @details assumes constant #IP/element
  */
void getSprPoints(Field* f,
                  std::set<MeshEntity*>& elements,
                  SprPoints& spr_points)
{
  Mesh* mesh = getMesh(f);
  std::set<MeshEntity*>::iterator elem = elements.begin();
  MeshElement* me = createMeshElement(mesh,*elem);
  int order = mesh->getShape()->getOrder();
  int num_points_elem = countIntPoints(me,order);
  destroyMeshElement(me);
  int np = elements.size()*num_points_elem;
  spr_points.setCount(np);
  std::size_t p = 0;
  for (std::set<MeshEntity*>::iterator it = elements.begin();
      it != elements.end(); ++it)
  {
    MeshElement* mesh_element = createMeshElement(mesh,*it);
    for (int l=0; l < num_points_elem; ++l)
    {
      Vector3 param;
      getIntPoint(mesh_element,order,l,param);
      mapLocalToGlobal(mesh_element,param,spr_points.points[p]);
      getMatrix(f,*it,l,spr_points.values[p]);
      ++p;
    }
    destroyMeshElement(mesh_element);
  }
}

/** @brief run patch recovery
  * @param f (In) integration point field
  * @param f_star (In) recovered nodal field
  * @elements (In) elements in patch
  * @param entity (In) central entity of patch
  */
void runSpr(Field* f,
            Field* f_star,
            std::set<MeshEntity*>& elements,
            MeshEntity* entity)
{
  SprPoints spr_points;
  getSprPoints(f,elements,spr_points);
  NewArray<double> component_values(spr_points.count);
  Vector3 point;
  getMesh(f)->getPoint(entity,0,point);
  Matrix3x3 value;
  Mesh* mesh = getMesh(f);
  int order = mesh->getShape()->getOrder();
  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
    {
      for (int p=0; p < spr_points.count; ++p)
        component_values[p] = spr_points.values[p][i][j];
      DynamicVector coeffs;
      evalPolynomialCoeffs(order,spr_points.count,
          spr_points.points,component_values,coeffs);
      DynamicVector terms;
      evalPolynomialTerms(order,point,terms);
      value[i][j] = coeffs*terms;
    }
  setMatrix(f_star,entity,0,value);
}

class PatchOp : public CavityOp
{
  public:
    PatchOp(Field* f_, Field* f_star_):
      CavityOp(getMesh(f_)),
      f(f_),
      f_star(f_star_),
      num_points(0)
    {
    }
    void reset(MeshEntity* e)
    {
      num_points = 0;
      elements.clear();
      entity = e;
    }
    virtual Outcome setEntity(MeshEntity* e)
    {
      if (hasEntity(f_star,e))
        return SKIP;
      reset(e);
      if ( ! buildPatch())
        return REQUEST;
      return OK;
    }
    virtual void apply()
    {
      runSpr(f,f_star,elements,entity);
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
      int o = mesh->getShape()->getOrder();
      int num_desired_points = (o+1)*(o+2)*(o+3)/6;
      if (num_points >= num_desired_points)
        return true;
      std::set<MeshEntity*> old_set = elements;
      int d = mesh->getDimension();
      for (int shared_dim=d-1; shared_dim >= 0; --shared_dim)
      {
        if (!addElementsThatShare(shared_dim,old_set))
          return false;
        if (num_points >= num_desired_points)
          return true;
      }
      bool hope = elements.size() > old_set.size();
      if (hope)
        return expandAsNecessary();
      else
      {
        fail("SPR: patch construction: all hope is lost.");
        return false;
      }
    }
    bool addElementsThatShare(int dimension,
                              std::set<MeshEntity*>& old_set)
    {
      std::set<MeshEntity*> bridges;
      for (std::set<MeshEntity*>::iterator it = old_set.begin();
          it != old_set.end(); ++it)
      {
        Downward down;
        int nd = mesh->getDownward(*it,dimension,down);
        for (int i=0; i < nd; ++i)
          bridges.insert(down[i]);
      }
      std::vector<MeshEntity*> bridge_array(bridges.begin(),bridges.end());
      bridges.clear();
      if ( ! requestLocality(&(bridge_array[0]),bridge_array.size()))
        return false;
      for (size_t i=0; i < bridge_array.size(); ++i)
      {
        Adjacent candidates;
        mesh->getAdjacent(bridge_array[i],mesh->getDimension(),candidates);
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
      int integration_order = f->getShape()->getOrder();
      num_points += countIntPoints(element,integration_order);
      destroyMeshElement(element);
    }
    Field* f;
    Field* f_star;
    int num_points;
    std::set<MeshEntity*> elements;
    MeshEntity* entity;
};

Field* recoverField(Field* f)
{
  std::string name = "spr_";
  name += getName(f);
  Mesh* m = getMesh(f);
  Field* f_star = createFieldOn(m,name.c_str(),MATRIX);
  PatchOp op(f,f_star);
  for (int i=0; i <= 3; ++i)
  {
    if (m->getShape()->countNodesOn(i) != 0)
      op.applyToDimension(i);
  }
  return f_star;
}

}

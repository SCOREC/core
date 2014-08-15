/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "spr.h"
#include "apfMesh.h"
#include "apfShape.h"
#include "apfField.h"
#include "apfCavityOp.h"
#include <set>

namespace spr {

struct SprPoints {
  SprPoints():num_points(0) {}
  void allocate(int np, int nc)
  {
    num_points = np;
    points.allocate(np);
    values.allocate(np);
    values2.allocate(np);
    for (int i=0; i < np; ++i)
      values2[i].allocate(nc);
  }
  int num_points;
  apf::NewArray<apf::Vector3> points;
  apf::NewArray<apf::Matrix3x3> values;
  apf::NewArray<apf::NewArray<double> > values2;
};

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
void evalPolynomialCoeffs(int order,
                          int num_points,
                          apf::NewArray<apf::Vector3> const& points,
                          apf::NewArray<double> const& values,
                          apf::DynamicVector& coeffs)
{
  int m = num_points;
  int n = (order+1)*(order+2)*(order+3)/6;
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

/** @brief get spr point data from element patch
  * @param f (In) integration point field
  * @param elements (In) elements in patch
  * @param spr_points (Out) coord and data at points
  * @details assumes constant #IP/element
  */
void getSprPoints(apf::Field* f,
                  std::set<apf::MeshEntity*>& elements,
                  SprPoints& spr_points)
{
  apf::Mesh* mesh = apf::getMesh(f);
  std::set<apf::MeshEntity*>::iterator elem = elements.begin();
  apf::MeshElement* me = apf::createMeshElement(mesh,*elem);
  int order = mesh->getShape()->getOrder();
  int num_points_elem = apf::countIntPoints(me,order);
  apf::destroyMeshElement(me);
  int np = elements.size()*num_points_elem;
  int nc = f->countComponents();
  spr_points.allocate(np,nc);
  std::size_t p = 0;
  for (std::set<apf::MeshEntity*>::iterator it = elements.begin();
      it != elements.end(); ++it)
  {
    apf::MeshElement* mesh_element = apf::createMeshElement(mesh,*it);
    for (int l=0; l < num_points_elem; ++l)
    {
      apf::Vector3 param;
      apf::getIntPoint(mesh_element,order,l,param);
      apf::mapLocalToGlobal(mesh_element,param,spr_points.points[p]);
      apf::getMatrix(f,*it,l,spr_points.values[p]);
      apf::getComponents(f,*it,l,&(spr_points.values2[p])[0]);
      ++p;
    }
    apf::destroyMeshElement(mesh_element);
  }
}

/** @brief run patch recovery
  * @param f (In) integration point field
  * @param f_star (In) recovered nodal field
  * @elements (In) elements in patch
  * @param entity (In) central entity of patch
  */
void runSpr(apf::Field* f,
            apf::Field* f_star,
            std::set<apf::MeshEntity*>& elements,
            apf::MeshEntity* entity)
{
  SprPoints spr_points;
  getSprPoints(f,elements,spr_points);
  apf::NewArray<double> component_values(spr_points.num_points);
  apf::Vector3 point;
  apf::getMesh(f)->getPoint(entity,0,point);
  apf::Matrix3x3 value;
  apf::Mesh* mesh = apf::getMesh(f);
  int order = mesh->getShape()->getOrder();
  for (int i=0; i < 3; ++i)
    for (int j=0; j < 3; ++j)
    {
      for (int p=0; p < spr_points.num_points; ++p)
        component_values[p] = spr_points.values[p][i][j];
      apf::DynamicVector coeffs;
      evalPolynomialCoeffs(order,spr_points.num_points,
          spr_points.points,component_values,coeffs);
      apf::DynamicVector terms;
      evalPolynomialTerms(order,point,terms);
      value[i][j] = coeffs*terms;
    }
  setMatrix(f_star,entity,0,value);
}

class PatchOp : public apf::CavityOp
{
  public:
    PatchOp(apf::Field* f_, apf::Field* f_star_):
      apf::CavityOp(apf::getMesh(f_)),
      f(f_),
      f_star(f_star_),
      num_points(0)
    {
    }
    void reset(apf::MeshEntity* e)
    {
      num_points = 0;
      elements.clear();
      entity = e;
    }
    virtual Outcome setEntity(apf::MeshEntity* e)
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
      apf::DynamicArray<apf::MeshEntity*> adjacent;
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
      std::set<apf::MeshEntity*> old_set = elements;
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
        apf::fail("SPR: patch construction: all hope is lost.");
        return false;
      }
    }
    bool addElementsThatShare(int dimension,
                              std::set<apf::MeshEntity*>& old_set)
    {
      std::set<apf::MeshEntity*> bridges;
      for (std::set<apf::MeshEntity*>::iterator it = old_set.begin();
          it != old_set.end(); ++it)
      {
        apf::Downward down;
        int nd = mesh->getDownward(*it,dimension,down);
        for (int i=0; i < nd; ++i)
          bridges.insert(down[i]);
      }
      std::vector<apf::MeshEntity*> bridge_array(bridges.begin(),bridges.end());
      bridges.clear();
      if ( ! requestLocality(&(bridge_array[0]),bridge_array.size()))
        return false;
      for (size_t i=0; i < bridge_array.size(); ++i)
      {
        apf::Adjacent candidates;
        mesh->getAdjacent(bridge_array[i],mesh->getDimension(),candidates);
        add(candidates);
      }
      return true;
    }
    void add(apf::DynamicArray<apf::MeshEntity*>& es)
    {
      for (std::size_t i=0; i < es.getSize(); ++i)
        add(es[i]);
    }
    void add(apf::MeshEntity* e)
    {
      if (elements.count(e))
        return;
      elements.insert(e);
      apf::MeshElement* element = apf::createMeshElement(mesh,e);
      int integration_order = f->getShape()->getOrder();
      num_points += apf::countIntPoints(element,integration_order);
      apf::destroyMeshElement(element);
    }
    apf::Field* f;
    apf::Field* f_star;
    int num_points;
    std::set<apf::MeshEntity*> elements;
    apf::MeshEntity* entity;
};

apf::Field* recoverField(apf::Field* f)
{
  std::string name = "spr_";
  name += getName(f);
  apf::Mesh* m = apf::getMesh(f);
  apf::Field* f_star = apf::createFieldOn(m,name.c_str(),apf::MATRIX);
  PatchOp op(f,f_star);
  for (int i=0; i <= 3; ++i)
  {
    if (m->getShape()->countNodesOn(i) != 0)
      op.applyToDimension(i);
  }
  return f_star;
}

}

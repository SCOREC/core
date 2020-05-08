/*
1;3409;0c * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include <apfElement.h>

#include "crv.h"
#include "crvAdapt.h"
#include "crvBezier.h"
#include "crvBezierShapes.h"
#include "crvMath.h"
#include "crvQuality.h"
#include "crvShape.h"
#include "crvTables.h"
#include <maMesh.h>
#include <maShapeHandler.h>
#include <maSolutionTransfer.h>
#include <maShape.h>
#include <maMap.h>
#include <mth_def.h>
#include <math.h>
#include <pcu_util.h>

#include <iostream>
namespace crv {

static void getVertParams(apf::Mesh* mesh, int ptype,
    apf::MeshEntity** parentVerts,
    apf::NewArray<apf::MeshEntity*>& midEdgeVerts,
    apf::MeshEntity* e,
    apf::Vector3 params[4])
{
  int npv = apf::Mesh::adjacentCount[ptype][0];
  int ne = apf::Mesh::adjacentCount[ptype][1];
  apf::Downward verts;

  int nv = mesh->getDownward(e,0,verts);
  // first check verts
  for (int v = 0; v < nv; ++v){
    bool vert = false;
    for (int i = 0; i < npv; ++i){
      if(verts[v] == parentVerts[i]){
        params[v] = elem_vert_xi[ptype][i];
        vert = true;
        break;
      }
    }

    if(!vert){
      for (int i = 0; i < ne; ++i){
        if( verts[v] == midEdgeVerts[i] ){
          params[v] = elem_edge_xi[ptype][i];
          break;
        }
      }
    }
  }
}

static int getBestElement(int n,
    apf::Mesh* mesh,
    apf::MeshElement **elems,
    apf::Vector3 point, apf::Vector3 &xi)
{
  int iter = 50;
  double tol = 1e-4;

  int elemNum = 0;
  double value, value2, bestValue = -DBL_MAX;

  apf::Vector3 initialGuess = apf::Vector3(1./4., 1./4., 1./4.);
  apf::Vector3 xinew, xyz;
  apf::Matrix3x3 Jinv;

  for (int i = 0; i < n; i++) {
    apf::Vector3 xin = initialGuess;

    for (int j = 0; j < iter; j++) {
      apf::getJacobianInv(elems[i], xin, Jinv);
      mapLocalToGlobal(elems[i], xin, xyz);
      apf::Vector3 fval = (xyz-point);

      xinew = xin - Jinv*fval;

      if (j > 0 && fval.getLength() > value) break;

      if ( (xinew-xin).getLength() < tol ) {
        value2 = ma::getInsideness(mesh, apf::getMeshEntity(elems[i]), xinew);
        if ( value2 > bestValue) {
          bestValue = value2;
          elemNum = i;
          xi = xinew;
          break;
        }
      }
      else {
        value = fval.getLength();
        xin = xinew;
      }
    }
  }
  return (elemNum);
}


class CrvBezierSolutionTransfer : public ma::SolutionTransfer
{
  public:
    Adapt* adapt;
    ma::LinearTransfer verts;
    ma::CavityTransfer others;
    CrvBezierSolutionTransfer(apf::Field* field, Adapt* a):
      adapt(a),verts(field),others(field)
    {
      f = field;
      refine = adapt->refine;
      mesh = adapt->mesh;
      shape = apf::getShape(f);
      int P = mesh->getShape()->getOrder();
      for (int d = 1; d <= mesh->getDimension(); d++) {
      	int type = apf::Mesh::simplexTypes[d];
      	if (!mesh->getShape()->hasNodesIn(d))
      	  continue;
      	int n = mesh->getShape()->getEntityShape(type)->countNodes();
	mth::Matrix<double> A(n,n);
	Ai[d].resize(n,n);
	getBezierTransformationMatrix(type, P,
	    A, elem_vert_xi[type]);
	invertMatrixWithPLU(getNumControlPoints(type,P), A, Ai[d]);
      }
    }
    virtual bool hasNodesOn(int dimension)
    {
      return others.hasNodesOn(dimension);
    }

    virtual void onVertex(
        apf::MeshElement* parent,
        ma::Vector const& xi,
        ma::Entity* vert)
    {
      verts.onVertex(parent,xi,vert);
    }

    virtual void onRefine(
        ma::Entity* parent,
        ma::EntityArray& newEntities)
    {
      int P = shape->getOrder();
      int parentType = mesh->getType(parent);
      apf::Downward parentVerts, parentEdges;

      mesh->getDownward(parent, 0, parentVerts);
      mesh->getDownward(parent, 1, parentEdges);

      int ne = apf::Mesh::adjacentCount[parentType][1];

      apf::NewArray<apf::MeshEntity*> midEdgeVerts(ne);
      for (int i = 0; i < ne; ++i){
        if ( ma::getFlag(adapt,parentEdges[i],ma::SPLIT) ) {
          midEdgeVerts[i] = ma::findSplitVert(refine,parentEdges[i]);
	  apf::Element* elemP = apf::createElement(f, parentEdges[i]);
	  apf::Vector3 xiMid = apf::Vector3(0,0,0);
	  apf::NewArray<double> val;
	  val.allocate(apf::countComponents(f));
	  apf::getComponents(elemP, xiMid, &(val[0]));
	  apf::setComponents(f, midEdgeVerts[i], 0, &(val[0]));
	  apf::destroyElement(elemP);
	}
        else
          midEdgeVerts[i] = 0;
      }

      int np = shape->getEntityShape(parentType)->countNodes();

      apf::Element* elem = apf::createElement(f, parent);
      apf::NewArray<apf::Vector3> Vnodes;
      apf::NewArray<double> Snodes;

      if (apf::getValueType(f) == apf::VECTOR) {
        apf::getVectorNodes(elem,Vnodes);
      }
      else if (apf::getValueType(f) == apf::SCALAR) {
	apf::getScalarNodes(elem,Snodes);
      }

      for (int d = 1; d <= apf::Mesh::typeDimension[parentType]; d++) {
        if (!shape->hasNodesIn(d)) continue;
        for (size_t i = 0; i < newEntities.getSize(); i++) {
          int childType = mesh->getType(newEntities[i]);
          if (d != apf::Mesh::typeDimension[childType])
            continue;

          int ni = shape->countNodesOn(childType);
          int n = getNumControlPoints(childType,P);
          apf::Vector3 vp[4];
          getVertParams(mesh, parentType,parentVerts,
              midEdgeVerts,newEntities[i], vp);

          mth::Matrix<double> A(n,np),B(n,n);
          getBezierTransformationMatrix(parentType, childType, P,A, vp);
          mth::multiply(Ai[apf::Mesh::typeDimension[childType]],A,B);

          for (int j = 0; j < ni; ++j){

            if (apf::getValueType(f) == apf::VECTOR) {
              apf::Vector3 Vvalue(0,0,0);
              for (int k = 0; k < np; ++k)
              	Vvalue += Vnodes[k]*B(j+n-ni,k);
              apf::setVector(f,newEntities[i],j,Vvalue);
	    }
	    else if (apf::getValueType(f) == apf::SCALAR) {
	      double Svalue = 0.;

	      for (int k = 0; k < np; ++k)
	      	Svalue += Snodes[k]*B(j+n-ni,k);
	      apf::setScalar(f,newEntities[i],j,Svalue);
	    }
          }
        }
      }
      apf::destroyElement(elem);
    }

    void setInterpolatingFieldValues(apf::Mesh* mesh,
    	ma::EntityArray &oldElements,
    	apf::MeshEntity* newEnt, int nodeNum,
    	apf::Vector3 xyz, apf::Field* field)
    {
      apf::Vector3 xi;
      int entNum = -1;
      apf::NewArray<apf::MeshElement*> elems(oldElements.getSize());
      for (size_t i = 0; i < oldElements.getSize(); i++)
      	elems[i] = apf::createMeshElement(mesh, oldElements[i]);
      // For a point outside the cavity getBestElement routine
      // returns an extrapolated xi
      // should perform a check for negative xi values
      entNum = getBestElement(oldElements.getSize(), mesh,
      	  &elems[0], xyz, xi);

      PCU_ALWAYS_ASSERT(entNum >= 0);

      apf::Element* elemP = apf::createElement(f, oldElements[entNum]);
      apf::NewArray<double> val;
      val.allocate(apf::countComponents(field));
      apf::getComponents(elemP, xi, &(val[0]));
      apf::setComponents(f, newEnt, nodeNum, &(val[0]));

      for (size_t i = 0; i < oldElements.getSize(); i++)
      	apf::destroyMeshElement(elems[i]);

      apf::destroyElement(elemP);
    }

    virtual void onCavity(
        ma::EntityArray& oldElements,
        ma::EntityArray& newEntities)
    {

      if (getDimension(mesh, oldElements[0]) < 3)
      	return;

      for (int d = 1; d < 4; d++) {
      	for(size_t i = 0; i <newEntities.getSize(); i++) {
      	  int type = mesh->getType(newEntities[i]);
      	  int td = apf::Mesh::typeDimension[type];
      	  if (td != d) continue;

      	  int order = shape->getOrder();
      	  int ne = shape->countNodesOn(apf::Mesh::simplexTypes[td]);

      	  if (!ne) continue;

      	  apf::MeshElement* me = apf::createMeshElement(mesh, newEntities[i]);
      	  apf::Vector3 point, xinew;
      	  int n = shape->getEntityShape(type)->countNodes();

      	  for (int j = 0; j < ne; j++) {
      	    shape->getNodeXi(type, j, xinew);
      	    apf::mapLocalToGlobal(me, xinew, point);
      	    setInterpolatingFieldValues(mesh, oldElements,
      	    	newEntities[i], j, point, f);
	  }
	  apf::NewArray<double> c;
	  crv::getBezierTransformationCoefficients(order,
	      apf::Mesh::simplexTypes[td], c);

	  crv::convertInterpolationFieldPoints(newEntities[i], f, n, ne, c);
	  apf::destroyMeshElement(me);
	}
      }
    }

  protected:
    apf::Field* f;
    ma::Mesh* mesh;
    ma::Refine* refine;
    apf::FieldShape* shape;
    mth::Matrix<double> Ai[4];
};

static ma::SolutionTransfer* createBezierSolutionTransfer(apf::Field* f, Adapt* a)
{
  return new CrvBezierSolutionTransfer(f, a);
}

ma::SolutionTransfer* setBezierSolutionTransfers(
    const std::vector<apf::Field*>& fields, Adapt* a)
{
  ma::SolutionTransfers* st = new ma::SolutionTransfers();
  for (std::size_t i = 0; i < fields.size(); i++) {
    st->add(createBezierSolutionTransfer(fields[i], a));
  }
  return st;
}

}

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
#include <string.h>
namespace crv {

static void setLinearEdgeFieldPoints(ma::Mesh* m, apf::Field* f,
    ma::Entity* edge)
{
  apf::Vector3 xi,points[2];
  apf::MeshEntity* verts[2];
  apf::FieldShape* fs = apf::getShape(f);

  int ni = fs->countNodesOn(apf::Mesh::EDGE);
  m->getDownward(edge,0,verts);
  apf::NewArray<double> value1, value2, value;
  value.allocate(apf::countComponents(f));
  value1.allocate(apf::countComponents(f));
  value2.allocate(apf::countComponents(f));
  apf::getComponents(f, verts[0], 0, &(value1[0]));
  apf::getComponents(f, verts[1], 0, &(value2[0]));

  m->getPoint(verts[0],0,points[0]);
  m->getPoint(verts[1],0,points[1]);
  for (int j = 0; j < ni; ++j){
    double t = (1.+j)/(1.+ni);
    for (int k = 0; k < apf::countComponents(f); k++)
      value[k] = value1[k]*(1.-t) + value2[k]*t;
    apf::setComponents(f, edge, j, &(value[0]));
  }
}

static void repositionInteriorFieldWithBlended(ma::Mesh* m,
    apf::Field* f, ma::Entity* e)
{
  apf::FieldShape * fs = apf::getShape(f);
  int order = fs->getOrder();
  int typeDim = apf::Mesh::typeDimension[m->getType(e)];
  int type = apf::Mesh::simplexTypes[typeDim];

  if(!fs->hasNodesIn(typeDim) || getBlendingOrder(type))
    return;

  int n = fs->getEntityShape(type)->countNodes();
  int ne = fs->countNodesOn(type);
  apf::NewArray<double> c;
  crv::getInternalBezierTransformationCoefficients(m,order,1,type,c);
  crv::convertInterpolationFieldPoints(e,f,n-ne,ne,c);

}


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
  int iter = 1000;
  double tol = 1e-4;

  int elemNum = 0;
  double value2, bestValue = -DBL_MAX;

  apf::Vector3 initialGuess = apf::Vector3(1./4., 1./4., 1./4.);
  apf::Vector3 xinew, xyz;
  apf::Matrix3x3 Jinv;
  apf::Vector3 fval;

  for (int i = 0; i < n; i++) {
    apf::Vector3 xin = initialGuess;
    apf::Element* e = apf::createElement(mesh->getCoordinateField(), elems[i]);

    for (int j = 0; j < iter; j++) {
      apf::getJacobianInv(elems[i], xin, Jinv);
      apf::getVector(e, xin, xyz);
      fval = (xyz-point);

      xinew = xin - transpose(Jinv)*fval;

      if ( (xinew-xin).getLength() < tol ) {
        value2 = ma::getInsideness(mesh, apf::getMeshEntity(elems[i]), xinew);
        if ( value2 > bestValue) {
          bestValue = value2;
          elemNum = i;
          xi = xinew;
          break;
        }
        else xin = xinew;
      }
      else {
        xin = xinew;
      }
    }
    apf::destroyElement(e);
  }

  return (elemNum);
}

class CrvBezierSolutionTransfer : public ma::SolutionTransfer
{
  public:
    Adapt* adapt;
    ma::CavityTransfer others;

    CrvBezierSolutionTransfer(apf::Field* field, Adapt* a):
      adapt(a),others(field)
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

    virtual const char* getTransferFieldName() {
      return apf::getName(f);
    }

    virtual void onVertex(
	apf::MeshElement* parent,
	ma::Vector const& xi,
	ma::Entity* vert)
    {
      apf::Element* e = apf::createElement(f, parent);
      apf::NewArray<double> value;
      value.allocate(apf::countComponents(f));
      apf::getComponents(e, xi, &(value[0]));
      apf::setComponents(f, vert, 0, &(value[0]));
      apf::destroyElement(e);
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
      bool useLinear = false;

      for (int i = 0; i < ne; ++i){
	if ( ma::getFlag(adapt,parentEdges[i],ma::SPLIT) )
	  midEdgeVerts[i] = ma::findSplitVert(refine,parentEdges[i]);
	else
	  midEdgeVerts[i] = 0;
	if ( ma::getFlag(adapt,parentEdges[i], ma::BAD_QUALITY) ) {
	  useLinear = true;
	}
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
	  int n = shape->getEntityShape(childType)->countNodes();
	  int ni = shape->countNodesOn(childType);

	  if (useLinear && !isBoundaryEntity(mesh, newEntities[i])) {
	    if (childType == apf::Mesh::EDGE) {
	      setLinearEdgeFieldPoints(mesh, f, newEntities[i]);
	    } else {
	      for (int j = 0; j < ni; j++) {
		apf::NewArray<double> val;
		val.allocate(apf::countComponents(f));
		for (int k = 0; k < apf::countComponents(f); k++)
		  val[k] = 0.;
		apf::setComponents(f, newEntities[i], j, &(val[0]));
	      }
	      repositionInteriorFieldWithBlended(mesh,f,newEntities[i]);
	    }
	  }
	  else {

	    apf::Vector3 vp[4];
	    getVertParams(mesh, parentType,parentVerts,
		midEdgeVerts,newEntities[i], vp);

	    mth::Matrix<double> A(n,np),B(n,np);
	    getBezierTransformationMatrix(parentType, childType, P,A, vp);
	    mth::multiply(Ai[apf::Mesh::typeDimension[childType]],A,B);

	    for (int j = 0; j < ni; ++j){

	      if (apf::getValueType(f) == apf::VECTOR) {
		apf::Vector3 Vvalue(0,0,0);
		for (int k = 0; k < np; ++k) {
		  Vvalue += Vnodes[k]*B(j+n-ni,k);
		}
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
      }
    }
    void setInterpolatingFieldValues(apf::Mesh* mesh,
	ma::EntityArray &oldElements,
	apf::MeshEntity* newEnt, int nodeNum,
	apf::Vector3 xyz, apf::Field* field1)
    {
      apf::Vector3 xi;
      int entNum = -1;
      int n = oldElements.getSize();

      apf::NewArray<apf::MeshElement*> elems(n);
      apf::NewArray<apf::Element*> elemsF(n);
      for (int i = 0; i < n; i++) {
	elems[i] = apf::createMeshElement(mesh, oldElements[i]);
      }

      entNum = getBestElement(n, mesh,
	  &elems[0], xyz, xi);

      PCU_ALWAYS_ASSERT(entNum >= 0);

      apf::MeshElement* meshElemP = apf::createMeshElement(
	  mesh, oldElements[entNum]);
      apf::Vector3 xxyy, xxzz, val2;
      apf::mapLocalToGlobal(meshElemP, xi, xxyy);

      //int parentType = mesh->getType(oldElements[entNum]);
      //int np = fshape->getEntityShape(parentType)->countNodes();

      apf::Element* elemP = apf::createElement(field1, oldElements[entNum]);
      apf::NewArray<double> val;
      val.allocate(apf::countComponents(field1));
      apf::getComponents(elemP, xi, &(val[0]));

      apf::setComponents(field1, newEnt, nodeNum, &(val[0]));
      for (int i = 0; i < n; i++)
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

	  int ne = shape->countNodesOn(apf::Mesh::simplexTypes[td]);

	  if (!ne) continue;

	  apf::MeshElement* me = apf::createMeshElement(mesh, newEntities[i]);
	  apf::Element* e = apf::createElement(mesh->getCoordinateField(), me);
	  apf::Vector3 point, xinew;

	  for (int j = 0; j < ne; j++) {
	    shape->getNodeXi(type, j, xinew);
	    apf::getVector(e, xinew, point);
	    //apf::mapLocalToGlobal(me, xinew, point);
	    setInterpolatingFieldValues(mesh, oldElements,
		newEntities[i], j, point, f);
	  }

	  apf::destroyMeshElement(me);
	}
      }
      // convert interpolating values to control points in the
      // decreasing order of entity dimension
      for (int d = 3; d >= 1; d--) {
	for(size_t i = 0; i < newEntities.getSize(); i++) {
	  int type = mesh->getType(newEntities[i]);
	  int td = apf::Mesh::typeDimension[type];
	  if (td != d) continue;

	  int order = shape->getOrder();
	  int ne = shape->countNodesOn(apf::Mesh::simplexTypes[td]);

	  if (!ne) continue;

	  int n = shape->getEntityShape(type)->countNodes();
	  apf::NewArray<double> c;
	  crv::getBezierTransformationCoefficients(order,
	      apf::Mesh::simplexTypes[td], c);

	  crv::convertInterpolationFieldPoints(newEntities[i], f, n, ne, c);
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
  ma::SolutionTransfers *st = dynamic_cast<ma::SolutionTransfers*>(a->solutionTransfer);
  if (!st)
    st = new ma::SolutionTransfers();

  for (std::size_t i = 0; i < fields.size(); i++) {
    st->add(createBezierSolutionTransfer(fields[i], a));
  }
  std::vector<ma::SolutionTransfer*> trans = st->transfers;
  for (std::size_t i = 0; i < trans.size(); i++) {
    const char* name = trans[i]->getTransferFieldName();
    std::cout<<" field name added "<< name<<std::endl;
  }
  return st;
}

}

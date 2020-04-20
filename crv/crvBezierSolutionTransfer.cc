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
#include <mth_def.h>
#include <math.h>
#include <pcu_util.h>

#include <iostream>
namespace crv {

class CrvBezierSolutionTransfer : public ma::SolutionTransfer
{
  public:
    Adapt* adapt;
    ma::LinearTransfer verts;
    ma::CavityTransfer others;
    CrvBezierSolutionTransfer(apf::Field* field, Adapt* a):
      adapt(a),verts(f),others(f)
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
	//crv::getBezierTransformationCoeddicients(P,type,coeffs[d]);
	//crv::getInternalBezierTransformationCoefficients(mesh, P, 1,
	//  apf::Mesh::simplexTypes[d],internalCoeffs[d]);
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

    void getVertParams(int ptype, apf::MeshEntity** parentVerts,
        apf::NewArray<apf::MeshEntity*>& midEdgeVerts, apf::MeshEntity* e,
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

    virtual void onRefine(
        ma::Entity* parent,
        ma::EntityArray& newEntities)
    {
      others.onRefine(parent,newEntities);

      int P = shape->getOrder();
      int parentType = mesh->getType(parent);
      apf::Downward parentVerts, parentEdges;

      mesh->getDownward(parent, 0, parentVerts);
      mesh->getDownward(parent, 1, parentEdges);

      int ne = apf::Mesh::adjacentCount[parentType][1];

      apf::NewArray<apf::MeshEntity*> midEdgeVerts(ne);
      for (int i = 0; i < ne; ++i){
        if ( ma::getFlag(adapt,parentEdges[i],ma::SPLIT) )
          midEdgeVerts[i] = ma::findSplitVert(refine,parentEdges[i]);
        else
          midEdgeVerts[i] = 0;
      }

      int np = shape->getEntityShape(parentType)->countNodes();

      apf::Element* elem = apf::createElement(f, parent);
      apf::NewArray<apf::Vector3> nodes;

      //TODO scalar field values
      apf::getVectorNodes(elem,nodes);

      for (int d = 1; d <= apf::Mesh::typeDimension[parentType]; d++) {
        if (!shape->hasNodesIn(d)) continue;
        for (size_t i = 0; i < newEntities.getSize(); i++) {
          int childType = mesh->getType(newEntities[i]);
          if (d != apf::Mesh::typeDimension[childType])
            continue;

          int ni = shape->countNodesOn(childType);
          int n = getNumControlPoints(childType,P);
          apf::Vector3 vp[4];
          getVertParams(parentType,parentVerts,
              midEdgeVerts,newEntities[i], vp);

          mth::Matrix<double> A(n,np),B(n,n);
          getBezierTransformationMatrix(parentType, childType, P,A, vp);
          mth::multiply(Ai[apf::Mesh::typeDimension[childType]],A,B);

          for (int j = 0; j < ni; ++j){
            apf::Vector3 value(0,0,0);
            for (int k = 0; k < np; ++k)
              value += nodes[k]*B(j+n-ni,k);
	    apf::setVector(f,newEntities[i],j,value);
          }
        }
      }
    }

    virtual void onCavity(
        ma::EntityArray& oldElements,
        ma::EntityArray& newEntities)
    {
      others.onCavity(oldElements,newEntities);
      apf::Mesh *m = others.mesh;
      for (size_t i = 0; i <newEntities.getSize(); i++) {
      	int type = m->getType(newEntities[i]);
      	apf::FieldShape *shape = others.shape;
      	int order = shape->getOrder();
      	if (type != 0 && getNumInternalControlPoints(type, order) > 0) {
      	  int td = apf::Mesh::typeDimension[type];
      	  int n = shape->getEntityShape(
      	      apf::Mesh::simplexTypes[td])->countNodes();
      	  int ne = shape->countNodesOn(
      	      apf::Mesh::simplexTypes[td]);
      	  apf::NewArray<double> c;
      	  crv::getBezierTransformationCoefficients(order,
      	      apf::Mesh::simplexTypes[td], c);

      	  crv::convertInterpolationFieldPoints(newEntities[i],
              others.field, n, ne, c);
        }
      }
    }
  protected:
    apf::Field* f;
    ma::Mesh* mesh;
    ma::Refine* refine;
    apf::FieldShape* shape;
    mth::Matrix<double> Ai[4];
  //public:
    //apf::NewArray<double> coeffs[4];
    //apf::NewArray<double> internalCoeffs[4];
    /*
    */
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

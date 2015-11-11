/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include <apfElement.h>

#include "crvAdapt.h"
#include "crvBezier.h"
#include "crvBezierShapes.h"
#include "crvMath.h"
#include "crvQuality.h"
#include "crvTables.h"
#include <maMap.h>
#include <maShapeHandler.h>
#include <maSolutionTransfer.h>
#include <maShape.h>
#include <mth_def.h>
#include <cassert>

namespace crv {

class BezierTransfer : public ma::SolutionTransfer
{
  public:
    BezierTransfer(ma::Mesh* m, ma::Refine* r, bool snap)
    {
      mesh = m;
      refine = r;
      shouldSnap = snap;
      // pre compute the inverses of the transformation matrices
      int P = mesh->getShape()->getOrder();
      for (int d = 1; d <= 3; ++d){
        if(!getNumInternalControlPoints(apf::Mesh::simplexTypes[d],P))
          continue;
        int n = getNumControlPoints(apf::Mesh::simplexTypes[d],P);
        mth::Matrix<double> A(n,n);
        Ai[d].resize(n,n);
        getBezierTransformationMatrix(apf::Mesh::simplexTypes[d],P,
            A,elem_vert_xi[apf::Mesh::simplexTypes[d]]);
        invertMatrixWithPLU(getNumControlPoints(apf::Mesh::simplexTypes[d],P),
            A,Ai[d]);
      }
    }
    ~BezierTransfer()
    {
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

        // this part relies on "closeness"
        // to determine if this is the correct edge
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
    virtual bool hasNodesOn(int dimension)
    {
      return mesh->getShape()->hasNodesIn(dimension);
    }
    virtual void onRefine(
        ma::Entity* parent,
        ma::EntityArray& newEntities)
    {
      int P = mesh->getShape()->getOrder();
      int parentType = mesh->getType(parent);

      // for the parent, get its vertices and mid edge nodes first
      apf::Downward parentVerts,parentEdges;

      mesh->getDownward(parent,0,parentVerts);
      mesh->getDownward(parent,1,parentEdges);
      int ne = apf::Mesh::adjacentCount[parentType][1];

      apf::NewArray<apf::MeshEntity*> midEdgeVerts(ne);
      for (int i = 0; i < ne; ++i)
        midEdgeVerts[i] = ma::findSplitVert(refine,parentEdges[i]);

      int np = getNumControlPoints(parentType,P);

      apf::Element* elem =
          apf::createElement(mesh->getCoordinateField(),parent);
      apf::NewArray<apf::Vector3> nodes;
      apf::getVectorNodes(elem,nodes);
      apf::destroyElement(elem);

      for (size_t i = 0; i < newEntities.getSize(); ++i)
      {
        int childType = mesh->getType(newEntities[i]);
        int ni = mesh->getShape()->countNodesOn(childType);

        if (childType == apf::Mesh::VERTEX || ni == 0 ||
            (mesh->getModelType(mesh->toModel(newEntities[i]))
            < mesh->getDimension() && shouldSnap))
          continue; //vertices will have been handled specially beforehand

        int n = getNumControlPoints(childType,P);

//        apf::Vector3 vp[4];
//        apf::Downward verts;
//        int nv = mesh->getDownward(newEntities[i],0,verts);
//        for (int j = 0; j < ni; ++j){
//          mesh->setPoint(newEntities[i],j,apf::Vector3(0,0,0));
//        }
//        apf::NewArray<apf::Vector3> oldNodes,newNodes(ni);
//        apf::Element* newElem =
//            apf::createElement(mesh->getCoordinateField(),newEntities[i]);
//        apf::getVectorNodes(newElem,oldNodes);
//        apf::destroyElement(newElem);
//
//        for (int v = 0; v < nv; ++v){
//          mesh->getPoint(verts[v],0,vp[v]);
//        }
//        for (int j = 0; j < ni; ++j){
//          apf::Vector3 xi;
//          mesh->getShape()->getNodeXi(childType,j,xi);
//          if(childType == apf::Mesh::EDGE){
//            xi[0] = 0.5*(xi[0]+1.);
//            oldNodes[j+n-ni] = vp[0]*(1.-xi[0])+vp[1]*xi[0];
//          }
//          if(childType == apf::Mesh::TRIANGLE){
//            oldNodes[j+n-ni] = vp[0]*(1.-xi[0]-xi[1])
//                + vp[1]*xi[0] + vp[2]*xi[1];
//          }
//          if(childType == apf::Mesh::TET){
//            oldNodes[j+n-ni] = vp[0]*(1.-xi[0]-xi[1]-xi[2]) + vp[1]*xi[0]
//                + vp[2]*xi[1] + vp[3]*xi[2];
//          }
//          mesh->setPoint(newEntities[i],j,oldNodes[j+n-ni]);
//        }
//        apf::NewArray<double> c;
//        crv::getBezierTransformationCoefficients(P,childType,c);
//        convertInterpolationPoints(n,ni,oldNodes,c,newNodes);
//        for (int j = 0; j < ni; ++j){
//          mesh->setPoint(newEntities[i],j,newNodes[j]);
//        }
        apf::Vector3 vp[4];
        getVertParams(parentType,parentVerts,midEdgeVerts,newEntities[i],vp);

        mth::Matrix<double> A(n,np),B(n,n);
        getBezierTransformationMatrix(parentType,childType,P,A,vp);
        mth::multiply(Ai[apf::Mesh::typeDimension[childType]],A,B);

        for (int j = 0; j < ni; ++j){
          apf::Vector3 point(0,0,0);
          for (int k = 0; k < np; ++k)
            point += nodes[k]*B(j+n-ni,k);

          mesh->setPoint(newEntities[i],j,point);
        }
      }
    }
  private:
    ma::Mesh* mesh;
    ma::Refine* refine;
    mth::Matrix<double> Ai[4];
    bool shouldSnap;
};

class BezierHandler : public ma::ShapeHandler
{
  public:
    BezierHandler(ma::Adapt* a)
    {
      mesh = a->mesh;
      bt = new BezierTransfer(mesh,a->refine,a->input->shouldSnap);
      sizeField = a->sizeField;
      shouldSnap = a->input->shouldSnap;
    }
    ~BezierHandler()
    {
      delete bt;
    }
    virtual double getQuality(apf::MeshEntity* e)
    {
      assert( mesh->getType(e) == apf::Mesh::TRIANGLE ||
          mesh->getType(e) == apf::Mesh::TET);
      return ma::measureElementQuality(mesh, sizeField, e)
        *crv::getQuality(mesh,e);
    }
    virtual bool hasNodesOn(int dimension)
    {
      return bt->hasNodesOn(dimension);
    }
    virtual void onRefine(
        apf::MeshEntity* parent,
        ma::EntityArray& newEntities)
    {
      bt->onRefine(parent,newEntities);
    }
    virtual void onCavity(
        ma::EntityArray& /*oldElements*/,
        ma::EntityArray& newEntities)
    {
      apf::FieldShape* fs = mesh->getShape();
      int P = fs->getOrder();
//      for (size_t j = 0; j < newEntities.getSize(); ++j){
//       apf::MeshEntity* edges[4];
//       int iEdge = 0;
//       apf::MeshEntity* newVerts[2];
//       mesh->getDownward(newEntities[j],0,newVerts);
//      // find the four edges that form the original hull
//      // each edge must have one vertex
//      for (size_t i = 0; i < oldElements.getSize(); ++i)
//      {
//        int oldType = mesh->getType(oldElements[i]);
//        if(oldType == apf::Mesh::EDGE){
//          apf::MeshEntity* verts[2];
//          mesh->getDownward(oldElements[i],0,verts);
//          if(newVerts[0] == verts[0] || newVerts[0] == verts[1]
//          || newVerts[0] == verts[0] || newVerts[0] == verts[1]){
//
//            iEdge++;
//          }
//        }
//      }
//
//    }

      // deal with all the boundary points, if a boundary edge has been
      // collapsed, this is a snapping operation
      int n = fs->getEntityShape(apf::Mesh::EDGE)->countNodes();
      apf::NewArray<double> c;
      crv::getBezierTransformationCoefficients(P,1,c);
      for (size_t i = 0; i < newEntities.getSize(); ++i)
      {
        int newType = mesh->getType(newEntities[i]);
        int ni = mesh->getShape()->countNodesOn(newType);
        if (newType != apf::Mesh::EDGE) continue;

        if (mesh->getModelType(mesh->toModel(newEntities[i]))
            < mesh->getDimension() && ni > 0 && shouldSnap){
            snapToInterpolate(mesh,newEntities[i]);
        } else {
          // use linearly spaced points
          apf::Vector3 xi,points[2];
          apf::MeshEntity* verts[2];
          mesh->getDownward(newEntities[i],0,verts);
          mesh->getPoint(verts[0],0,points[0]);
          mesh->getPoint(verts[1],0,points[1]);
          for (int j = 0; j < ni; ++j){
            fs->getNodeXi(apf::Mesh::EDGE,j,xi);
            double t = (xi[0]+1.)*0.5;
            xi = points[0]*(1.-t)+points[1]*t;
            mesh->setPoint(newEntities[i],j,xi);
          }
        }
        convertInterpolationPoints(mesh,newEntities[i],n,ni,c);
      }

      for (int d = 2; d <= mesh->getDimension(); ++d){

        int ni = fs->countNodesOn(apf::Mesh::simplexTypes[d]);
        if(ni == 0) continue;

        n = fs->getEntityShape(apf::Mesh::simplexTypes[d])->countNodes();
        apf::NewArray<double> c;
        crv::getInternalBezierTransformationCoefficients(mesh,
            P,1,apf::Mesh::simplexTypes[d],c);

        for (size_t i = 0; i < newEntities.getSize(); ++i)
        {
          int newType = mesh->getType(newEntities[i]);
          // zero the newEntities.
          if(apf::Mesh::typeDimension[newType] == d){

            for (int j = 0; j < ni; ++j){
              apf::Vector3 zero(0,0,0);
              mesh->setPoint(newEntities[i],j,zero);
            }
          }
          if (apf::Mesh::typeDimension[newType] == d && ni > 0
              && (mesh->getModelType(mesh->toModel(newEntities[i]))
              == mesh->getDimension() || !shouldSnap)){
            convertInterpolationPoints(mesh,newEntities[i],n-ni,ni,c);
          }
        }
      }
    }
  private:
    ma::Mesh* mesh;
    BezierTransfer* bt;
    ma::SizeField * sizeField;
    bool shouldSnap;
};

ma::ShapeHandler* getShapeHandler(ma::Adapt* a)
{
  return new BezierHandler(a);
}

}

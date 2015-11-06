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
    BezierTransfer(ma::Mesh* m, ma::Refine* r)
    {
      mesh = m;
      refine = r;
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

        if (childType == apf::Mesh::VERTEX || ni == 0)
          continue; //vertices will have been handled specially beforehand

        int n = getNumControlPoints(childType,P);

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
};

class BezierHandler : public ma::ShapeHandler
{
  public:
    BezierHandler(ma::Adapt* a)
    {
      mesh = a->mesh;
      bt = new BezierTransfer(mesh,a->refine);
      ct = ma::createFieldTransfer(mesh->getCoordinateField());
      sizeField = a->sizeField;
    }
    ~BezierHandler()
    {
      delete bt;
      delete ct;
    }
    virtual double getQuality(apf::MeshEntity* e)
    {
      assert( mesh->getType(e) == apf::Mesh::TRIANGLE ||
          mesh->getType(e) == apf::Mesh::TET);
      return crv::getQuality(mesh,e)*
          ma::measureElementQuality(mesh, sizeField, e);
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
        ma::EntityArray& oldElements,
        ma::EntityArray& newEntities)
    {
      ct->onCavity(oldElements,newEntities);
    }
  private:
    ma::Mesh* mesh;
    BezierTransfer* bt;
    ma::SolutionTransfer* ct;
    ma::SizeField * sizeField;

};

ma::ShapeHandler* getShapeHandler(ma::Adapt* a)
{
  return new BezierHandler(a);
}

}

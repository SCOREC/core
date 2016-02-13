/*
1;3409;0c * Copyright 2015 Scientific Computation Research Center
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
#include <maShapeHandler.h>
#include <maSolutionTransfer.h>
#include <maShape.h>
#include <mth_def.h>
#include <math.h>
#include <cassert>

namespace crv {

static double measureLinearTriArea(ma::Mesh* m, ma::Entity* tri)
{
  ma::Vector p[3];
  ma::getVertPoints(m,tri,p);
  return 0.5*apf::cross(p[1]-p[0],p[2]-p[0]).getLength();
}

static void setLinearEdgePoints(ma::Mesh* m, ma::Entity* edge)
{
  apf::Vector3 xi,points[2];
  apf::MeshEntity* verts[2];
  int ni = m->getShape()->countNodesOn(apf::Mesh::EDGE);
  m->getDownward(edge,0,verts);
  m->getPoint(verts[0],0,points[0]);
  m->getPoint(verts[1],0,points[1]);
  for (int j = 0; j < ni; ++j){
    double t = (1.+j)/(1.+ni);
    xi = points[0]*(1.-t)+points[1]*t;
    m->setPoint(edge,j,xi);
  }
}

class BezierTransfer : public ma::SolutionTransfer
{
  public:
    BezierTransfer(ma::Adapt* a)
    {
      adapt = a;
      mesh = a->mesh;
      refine = a->refine;
      shouldSnap = a->input->shouldSnap;
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

      // check if we can use the curvature of the original element or not
      // uses BAD_QUALITY on edges to indicate this
      bool useLinear = false;

      apf::NewArray<apf::MeshEntity*> midEdgeVerts(ne);
      for (int i = 0; i < ne; ++i){
        if ( ma::getFlag(adapt,parentEdges[i],ma::SPLIT) )
          midEdgeVerts[i] = ma::findSplitVert(refine,parentEdges[i]);
        else
          midEdgeVerts[i] = 0;
        if ( ma::getFlag(adapt,parentEdges[i],ma::BAD_QUALITY) ){
          useLinear = true;
        }
      }
      // do snapping inside here, since its not done elsewhere
      // general assumption is that initial mesh interpolates the geometry
      // well enough that we can snap and not affect the quality/validity
      // at a "significant level"
      // Since split vertex is already on the shape
      // of the parent, snapping should not move it very far
      //
      // do vertices first
      for (int i = 0; i < ne; ++i){
        if(shouldSnap && midEdgeVerts[i] &&
            isBoundaryEntity(mesh,midEdgeVerts[i])){
          snapToInterpolate(mesh,midEdgeVerts[i]);
        }
      }
      int np = getNumControlPoints(parentType,P);

      apf::Element* elem =
          apf::createElement(mesh->getCoordinateField(),parent);
      apf::NewArray<apf::Vector3> nodes;
      apf::getVectorNodes(elem,nodes);
      apf::destroyElement(elem);

      for (int d = 1; d <= apf::Mesh::typeDimension[parentType]; ++d){
        if (!mesh->getShape()->hasNodesIn(d)) continue;
        for (size_t i = 0; i < newEntities.getSize(); ++i)
        {
          // go through this hierachically, doing edges first
          int childType = mesh->getType(newEntities[i]);
          if(d != apf::Mesh::typeDimension[childType])
            continue;

          int ni = mesh->getShape()->countNodesOn(childType);
          bool isBdryEnt = isBoundaryEntity(mesh,newEntities[i]);

          // do snapping here, inside refinement
          if (isBdryEnt && shouldSnap){
            snapToInterpolate(mesh,newEntities[i]);
          } else if (useLinear && !isBdryEnt) {
            // boundary entities that don't get snapped should interpolate
            if(childType == apf::Mesh::EDGE){
              setLinearEdgePoints(mesh,newEntities[i]);
            } else {
              for (int j = 0; j < ni; ++j){
                apf::Vector3 zero(0,0,0);
                mesh->setPoint(newEntities[i],j,zero);
              }
              repositionInteriorWithBlended(mesh,newEntities[i]);
            }
          } else {
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
      }
      // again have to do this hierarchically, downward
      for (int d = apf::Mesh::typeDimension[parentType]; d >= 1; --d){
        if (!mesh->getShape()->hasNodesIn(d)) continue;
        for (size_t i = 0; i < newEntities.getSize(); ++i)
        {
          // go through this hierachically, doing edges first
          int childType = mesh->getType(newEntities[i]);
          if(d != apf::Mesh::typeDimension[childType] ||
              !isBoundaryEntity(mesh,newEntities[i]) || !shouldSnap)
            continue;

          int n = mesh->getShape()->getEntityShape(apf::Mesh::simplexTypes[d])
              ->countNodes();
          int ni = mesh->getShape()->countNodesOn(d);
          apf::NewArray<double> c;
          crv::getBezierTransformationCoefficients(P,d,c);
          convertInterpolationPoints(mesh,newEntities[i],n,ni,c);
        }
      }
    }
  private:
    ma::Adapt* adapt;
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
      adapt = a;
      mesh = a->mesh;
      bt = new BezierTransfer(a);
      sizeField = a->sizeField;
      shouldSnap = a->input->shouldSnap;
    }
    ~BezierHandler()
    {
      delete bt;
    }
    virtual double getQuality(apf::MeshEntity* e)
    {
      if (mesh->getType(e) == apf::Mesh::TRIANGLE){
        ma::Vector p[3];
        ma::getVertPoints(mesh,e,p);
        double l[3];
        for (int i=0; i < 3; ++i)
          l[i] = (p[(i+1)%3]-p[i]).getLength();
        double A = 0.5*apf::cross(p[1]-p[0],p[2]-p[0])[2];
        double s=0;
        for (int i=0; i < 3; ++i)
          s += l[i]*l[i];
        double lq;
        if (A < 0)
          lq = -48*(A*A)/(s*s);
        else
          lq = 48*(A*A)/(s*s);
        if (lq < 0)
          return lq;
        else return lq*crv::getQuality(mesh,e);
      }
      if (mesh->getType(e) == apf::Mesh::TET){
        ma::Vector p[4];
        ma::getVertPoints(mesh,e,p);
        double lq = ma::measureLinearTetQuality(p);
        if (lq < 0)
          return lq;
        else return lq*crv::getQuality(mesh,e);
      }
      return -1;
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

    // for a given edge and a cavity, attempt to find the pair
    // of triangles that the edge spans. This is not always possible
    // if multiple choices exist, as in 3D,
    // then use the pair with the minimum area
    bool findEdgeTrianglesCross(ma::EntityArray& entities,
        ma::Entity* edgeVerts[2],
        ma::Entity* edgeTris[2])
    {
      ma::Entity* edge0 = 0;
      ma::Entity* edge1 = 0;
      edgeTris[0] = edgeTris[1] = 0;
      double totalArea = 1e10;
      for (size_t i = 0; i < entities.getSize(); ++i){
        if (mesh->getType(entities[i]) != apf::Mesh::TRIANGLE) continue;
        // find a triangle with the first vertex of the edge
        if (ma::isInClosure(mesh,entities[i],edgeVerts[0])){
          // try the opposite edge, see if there's a match opposite it
          edge0 = ma::getTriEdgeOppositeVert(mesh,entities[i],edgeVerts[0]);
          for (size_t j = 0; j < entities.getSize(); ++j){
            if(j == i) continue;
            if(ma::isInClosure(mesh,entities[j],edgeVerts[1])){
              edge1 = ma::getTriEdgeOppositeVert(mesh,entities[j],edgeVerts[1]);
              if(edge0 == edge1){
                double newTotalArea = measureLinearTriArea(mesh, entities[i])
                  + measureLinearTriArea(mesh,entities[j]);
                if(newTotalArea < totalArea){
                  edgeTris[0] = entities[i];
                  edgeTris[1] = entities[j];
                  totalArea = newTotalArea;
                }
              }
            }
          }
        }
      }
      if(edgeTris[0])
        return true;
      return false;
    }
    // for a given edge and a cavity, attempt to find the pair
    // of triangles that share the edge
    // in 2D this is trivial,
    // in 3D it requires a choice
    // The edge we are trying to find triangles around
    // has four triangles to choose from
    // pick the two smallest based on linear area,
    // since we don't know any better
    bool findEdgeTrianglesShared(ma::Entity* edge,
        ma::Entity* edgeTris[2])
    {
      apf::Up up;
      mesh->getUp(edge,up);
      if(mesh->getDimension() == 2){
        edgeTris[0] = up.e[0];
        edgeTris[1] = up.e[1];
      } else {
        edgeTris[0] = edgeTris[1] = 0;
        double a0, a1;
        a0 = a1 = 1e10;
        ma::Entity* verts[2];
        mesh->getDownward(edge,0,verts);
        for (int i = 0; i < up.n; ++i){
          double a = measureLinearTriArea(mesh,up.e[i]);
          if(a < a0){
            a1 = a0;
            a0 = a;
            edgeTris[1] = edgeTris[0];
            edgeTris[0] = up.e[i];
          } else if (a < a1){
            a1 = a;
            edgeTris[1] = up.e[i];
          }
        }
      }
      assert(edgeTris[0] && edgeTris[1]);
      return true;
    }
    void evaluateBlendedQuad(ma::Entity* verts[4], ma::Entity* edges[4],
        int dir[4], apf::Vector3& xi, apf::Vector3& point)
    {
      point.zero();
      apf::Vector3 pt;
      apf::Vector3 xii[4] = {
          apf::Vector3(2.*xi[0]-1.,0,0),apf::Vector3(2.*xi[1]-1.,0,0),
          apf::Vector3(1.-2.*xi[0],0,0),apf::Vector3(1.-2.*xi[1],0,0)};
      // coefficients
      double eC[4] = {1.-xi[1],xi[0],xi[1],1.-xi[0]};
      double vC[4] = {eC[3]*eC[0],eC[0]*eC[1],eC[1]*eC[2],eC[2]*eC[3]};
      for (int i = 0; i < 4; ++i){
        apf::Element* elem =
            apf::createElement(mesh->getCoordinateField(),edges[i]);
        if(!dir[i]) xii[i][0] *= -1.;
        apf::getVector(elem,xii[i],pt);
        point += pt*eC[i]-ma::getPosition(mesh,verts[i])*vC[i];
        apf::destroyElement(elem);
      }
    }

    void setBlendedQuadEdgePoints(ma::Entity* edge,
        ma::Entity* verts[4], ma::Entity* edges[4],
        int dir[4])
    {
      int P = mesh->getShape()->getOrder();
      if(P == 2){
        apf::Vector3 xi(0.5,0.5,0);
        apf::Vector3 point;
        evaluateBlendedQuad(verts,edges,dir,xi,point);
        mesh->setPoint(edge,0,point);
      } else {
        apf::Vector3 xi(1./3.,1./3.,0);
        apf::Vector3 point;
        evaluateBlendedQuad(verts,edges,dir,xi,point);
        mesh->setPoint(edge,0,point);
        xi[0] = 2./3.; xi[1] = 2./3.;
        evaluateBlendedQuad(verts,edges,dir,xi,point);
        mesh->setPoint(edge,1,point);
        if (P > 3)
          elevateBezierCurve(mesh,edge,3,P-3);
      }
//      for (int i = 0; i < P-1; ++i)
//      {
//        double x = (1.+i)/P;
//        apf::Vector3 xi(x,x,0);
//        apf::Vector3 point;
//        evaluateBlendedQuad(verts,edges,dir,xi,point);
//        mesh->setPoint(edge,i,point);
//      }
    }
    ma::Entity* findEdgeInTri(ma::Entity* v0, ma::Entity* v1,
        ma::Entity* tri, int& dir)
    {
      ma::Entity* edges[3];
      mesh->getDownward(tri,1,edges);
      for (int i = 0; i < 3; ++i){
        ma::Entity* verts[2];
        mesh->getDownward(edges[i],0,verts);
        if(verts[0] == v0 && verts[1] == v1){
          dir = 1;
          return edges[i];
        }
        if(verts[0] == v1 && verts[1] == v0){
          dir = 0;
          return edges[i];
        }
      }
      fail("can't find edge in tri");
    }
    bool setBlendedQuadEdgePointsCross(ma::EntityArray& cavity,
        ma::Entity* edge)
    {
      ma::Entity* edgeVerts[2];
      ma::Entity* edgeTris[2];
      mesh->getDownward(edge,0,edgeVerts);
      if(findEdgeTrianglesCross(cavity,edgeVerts,edgeTris)){
        int index[4];
        ma::Entity* vA[3], *vB[3];
        mesh->getDownward(edgeTris[0],0,vA);
        mesh->getDownward(edgeTris[1],0,vB);
        index[0] = apf::findIn(vA,3,edgeVerts[0]);
        index[1] = apf::findIn(vB,3,vA[(index[0]+1) % 3]);
        index[2] = apf::findIn(vB,3,edgeVerts[1]);
        index[3] = (index[0]+2) % 3;

        ma::Entity* verts[4] = {vA[index[0]],vB[index[1]],
            vB[index[2]],vA[index[3]]};

        ma::Entity* edges[4];
        int dir[4];
        for (int i = 0; i < 4; ++i)
          edges[i] = findEdgeInTri(verts[i],verts[(i+1)%4],
              edgeTris[i==1||i==2],dir[i]);
        setBlendedQuadEdgePoints(edge,verts,edges,dir);
        return true;
      }
      return false;
    }
    void setBlendedQuadEdgePointsShared(ma::Entity* edge)
    {
      ma::Entity* edgeVerts[2];
      ma::Entity* edgeTris[2];
      mesh->getDownward(edge,0,edgeVerts);
      findEdgeTrianglesShared(edge,edgeTris);
      int index[4];
      ma::Entity* vA[3], *vB[3];
      mesh->getDownward(edgeTris[0],0,vA);
      mesh->getDownward(edgeTris[1],0,vB);
      index[0] = apf::findIn(vA,3,edgeVerts[0]);
      index[1] = apf::findIn(vA,3,
          ma::getTriVertOppositeEdge(mesh,edgeTris[0],edge));
      index[2] = apf::findIn(vB,3,edgeVerts[1]);
      index[3] = apf::findIn(vB,3,
          ma::getTriVertOppositeEdge(mesh,edgeTris[1],edge));
      assert(index[0] >= 0 && index[1] >= 0);
      assert(index[2] >= 0 && index[3] >= 0);
      ma::Entity* verts[4] = {vA[index[0]], vA[index[1]],
          vB[index[2]],vB[index[3]]};
      ma::Entity* edges[4];
      int dir[4];
      for (int i = 0; i < 4; ++i)
        edges[i] = findEdgeInTri(verts[i],verts[(i+1)%4],edgeTris[i>1],dir[i]);

      setBlendedQuadEdgePoints(edge,verts,edges,dir);
    }
    virtual void onCavity(
        ma::EntityArray& oldElements,
        ma::EntityArray& newEntities)
    {
      apf::FieldShape* fs = mesh->getShape();
      int P = fs->getOrder();
      int n = fs->getEntityShape(apf::Mesh::EDGE)->countNodes();
      apf::NewArray<double> edgeC;
      apf::NewArray<double> triC;

      crv::getBezierTransformationCoefficients(P,1,edgeC);
      if (P > 2)
        crv::getBezierTransformationCoefficients(P,2,triC);

      int numNewElements = 0;
      int numNewTriangles = 0;
      int numMiddleEdges = 0; // upper bound

      // deal with all the boundary points, if a boundary edge has been
      // collapsed, this is a snapping operation
      // also count a few things for later use, here
      ma::EntityArray oldTriangles;
      if(mesh->getDimension() == 2){
        oldTriangles = oldElements;
      } else if(mesh->getDimension() == 3){
        // find the set of triangles on the inside
        ma::EntitySet triangleSet;
        ma::Entity* tris[4];
        for (size_t i = 0; i < oldElements.getSize(); ++i)
        {
          mesh->getDownward(oldElements[i],2,tris);
          for (int j = 0; j < 4; ++j)
          {
            if(!triangleSet.insert(tris[j]).second)
              oldTriangles.append(tris[j]);
          }
        }
      }
      for (size_t i = 0; i < newEntities.getSize(); ++i)
      {
        int newType = mesh->getType(newEntities[i]);
        int ni = mesh->getShape()->countNodesOn(newType);

        if (newType == apf::Mesh::TRIANGLE)
          numNewTriangles++;
        if (newType == apf::Mesh::TET)
          numNewElements++;
        // zero new entities
        bool snap = isBoundaryEntity(mesh,newEntities[i]) && shouldSnap;
        if (newType != apf::Mesh::EDGE)
        {
          if (snap && P > 2){
            snapToInterpolate(mesh,newEntities[i]);
            convertInterpolationPoints(mesh,newEntities[i],n,ni,triC);
          } else {
            for (int j = 0; j < ni; ++j){
              apf::Vector3 zero(0,0,0);
              mesh->setPoint(newEntities[i],j,zero);
            }
          }
        } else {
          if (snap){
            snapToInterpolate(mesh,newEntities[i]);
            convertInterpolationPoints(mesh,newEntities[i],n,ni,edgeC);
          } else {
            setLinearEdgePoints(mesh,newEntities[i]);
            numMiddleEdges++;
          }
        }
      }

      ma::EntityArray middleEdges(numMiddleEdges);
      int me = 0;
      for (size_t i = 0; i < newEntities.getSize(); ++i){
        // if we aren't an edge or we are on a boundary, don't do this
        if(mesh->getType(newEntities[i]) != apf::Mesh::EDGE) continue;
        if(isBoundaryEntity(mesh,newEntities[i])) continue;
        // special case in 2D
        if(numNewTriangles == 2 && mesh->getDimension() == 2)
          setBlendedQuadEdgePointsShared(newEntities[i]);

        else if(!setBlendedQuadEdgePointsCross(oldTriangles,newEntities[i])){
          middleEdges[me] = newEntities[i];
          me++;
          // set to linear, because we don't know any better
          setLinearEdgePoints(mesh,newEntities[i]);
        }
      }

      // set the middle edges
      for (int i = 0; i < me; ++i)
        setBlendedQuadEdgePointsShared(middleEdges[i]);

      // set the rest of the interior points
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
          bool boundary = isBoundaryEntity(mesh,newEntities[i]);
          if (apf::Mesh::typeDimension[newType] == d && ni > 0
              && (!boundary || !shouldSnap)){
            convertInterpolationPoints(mesh,newEntities[i],n-ni,ni,c);
          }
        }
      }
    }
  private:
    ma::Adapt* adapt;
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

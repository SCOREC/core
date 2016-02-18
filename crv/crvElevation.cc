/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crv.h"
#include "crvBezier.h"
#include "crvMath.h"
#include "crvShape.h"
#include "crvSnap.h"
#include "crvTables.h"
#include "crvQuality.h"
#include <apfTagData.h>
#include <apfVectorField.h>

namespace crv {

void changeMeshOrder(apf::Mesh2* m, int newOrder)
{
  std::string name = m->getShape()->getName();
  if(name != std::string("Bezier"))
    fail("mesh must be already bezier");

  int oldOrder = m->getShape()->getOrder();
  if(oldOrder == newOrder)
   return;

  apf::VectorField* newCoordinateField = new apf::VectorField();
  newCoordinateField->init("__new_coordinates",
      m, crv::getBezier(newOrder), new apf::TagDataOf<double>());

  // do vertices first
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(0);
  apf::Vector3 coord;
  while ((e = m->iterate(it))) {
    m->getPoint(e,0,coord);
    apf::setVector(newCoordinateField,e,0,coord);
  }
  m->end(it);

  apf::Vector3 xi;

  // currently, we don't allow for variable order meshes,
  // as such, the bezier shape class can only be one order at a time
  // to make things happy, we switch back and forth depending on
  // whether we are elevating elements or snapping to boundary
  setOrder(oldOrder);
  bool canSnap = m->canSnap();

  // do the boundaries first for the new field, this is tricky,
  // snapping them

  for (int d = 1; d <= 3; ++d){
    it = m->begin(d);
    int type = apf::Mesh::simplexTypes[d];
    int nNewOn = getNumInternalControlPoints(type,newOrder);
    apf::Vector3 p;
    while ((e = m->iterate(it))) {
      if(canSnap && isBoundaryEntity(m,e)){
        for(int i = 0; i < nNewOn; ++i){
          apf::ModelEntity* g = m->toModel(e);
          getBezierNodeXi(type,newOrder,i,xi);
          if(type == apf::Mesh::EDGE)
            transferParametricOnEdgeSplit(m,e,0.5*(xi[0]+1.),p);
          else
            transferParametricOnTriSplit(m,e,xi,p);
          m->snapToModel(g,p,coord);
          apf::setVector(newCoordinateField,e,i,coord);
        }
      } else if (newOrder < oldOrder) {
        // decrease the order, using old mesh
        apf::Element* oldElem = apf::createElement(m->getCoordinateField(),e);
        for (int i = 0; i < nNewOn; ++i){
          getBezierNodeXi(type,newOrder,i,xi);
          apf::getVector(oldElem,xi,coord);
          apf::setVector(newCoordinateField,e,i,coord);
        }
        apf::destroyElement(oldElem);
      }
    }
    m->end(it);
  }
  setOrder(newOrder);

  // then go downward, and change the internal entities
  for (int d = m->getDimension(); d >= 1; --d){
    int type = apf::Mesh::simplexTypes[d];
    int nNewOn = getNumInternalControlPoints(type,newOrder);
    int nNew = getNumControlPoints(type,newOrder);
    if (!nNewOn) continue;

    apf::NewArray<apf::Vector3> oldNodes;
    apf::NewArray<apf::Vector3> newNodes(nNew);
    it = m->begin(d);

    // offset used in elevation
    int offset = nNew-nNewOn;
    if (type == apf::Mesh::EDGE)
      offset = 1;

    apf::NewArray<double> c;
    crv::getBezierTransformationCoefficients(newOrder,type,c);

    while ((e = m->iterate(it))) {
      if(newOrder < oldOrder || (canSnap && isBoundaryEntity(m,e))){
        // create element to change interpolating points
        // to control points, which are either on the boundary, or everywhere
        // depending on the case
        apf::Element* newElem = apf::createElement(newCoordinateField,e);
        apf::getVectorNodes(newElem,oldNodes);
        convertInterpolationPoints(nNew,nNewOn,oldNodes,c,newNodes);
        for(int i = 0; i < nNewOn; ++i)
          apf::setVector(newCoordinateField,e,i,newNodes[i]);
        apf::destroyElement(newElem);
      } else {
        // elevate the order to create elements from coordinate field
        setOrder(oldOrder);
        apf::Element* oldElem = apf::createElement(m->getCoordinateField(),e);
        apf::getVectorNodes(oldElem,oldNodes);
        elevateBezier(type,oldOrder,newOrder-oldOrder,oldNodes,newNodes);
        setOrder(newOrder); //set order back
        for(int i = 0; i < nNewOn; ++i){
          apf::setVector(newCoordinateField,e,i,newNodes[offset+i]);
        }
        apf::destroyElement(oldElem);

      }
    }
    m->end(it);
  }
  setOrder(newOrder);
  // set the coordinate field to the newly created one
  m->setCoordinateField(newCoordinateField);
}

/*
 * Templating is used for coordinates (Vector3) and det(Jacobian) (double)
 * and is only accessible in this file.
 */
template <class T>
static void raiseBezierEdge(int P, int r, apf::NewArray<T>& nodes,
    apf::NewArray<T>& elevatedNodes)
{
  elevatedNodes[0] = nodes[0];
  elevatedNodes[P+r] = nodes[P];
  for(int i = 1; i < P+r; ++i){
    for(int j = std::max(0,i-r); j <= std::min(i,P); ++j)
      elevatedNodes[i] += nodes[j]*binomial(P,j)*binomial(r,i-j)
      /binomial(P+r,i);
  }
}

void elevateBezierEdge(int P, int r, apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3>& elevatedNodes)
{
  // re-order nodes, makes life easier
  apf::Vector3 temp = nodes[1];
  for (int i = 1; i < P; ++i)
    nodes[i] = nodes[i+1];
  nodes[P] = temp;

  for (int i = 1; i < P+r; ++i)
    elevatedNodes[i].zero();
  raiseBezierEdge(P,r,nodes,elevatedNodes);
}

static void elevateBezierEdgeJacobianDet(int P, int r,
    apf::NewArray<double>& nodes,
    apf::NewArray<double>& elevatedNodes)
{
  for (int i = 1; i < P+r; ++i)
    elevatedNodes[i] = 0.;
  raiseBezierEdge(P,r,nodes,elevatedNodes);
}

template <class T>
static void raiseBezierTriangle(int P, int r, apf::NewArray<T>& nodes,
    apf::NewArray<T>& elevatedNodes)
{
  for(int i = 0; i <= P+r; ++i){
    for(int j = 0; j <= P+r-i; ++j){
      for(int k = std::max(0,i-r); k <= std::min(i,P); ++k){
        for(int l = std::max(0,i-k+j-r); l <= std::min(j,P-k); ++l){
          elevatedNodes[computeTriNodeIndex(P+r,i,j)] +=
              nodes[computeTriNodeIndex(P,k,l)]*trinomial(P,k,l)
              *trinomial(r,i-k,j-l)/trinomial(P+r,i,j);
        }
      }
    }
  }
}

void elevateBezierTriangle(int P, int r, apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3>& elevatedNodes)
{
  for(int i = 0; i < getNumControlPoints(apf::Mesh::TRIANGLE,P+r); ++i)
    elevatedNodes[i].zero();
  raiseBezierTriangle(P,r,nodes,elevatedNodes);
}

static void elevateBezierTriangleJacobianDet(int P, int r,
    apf::NewArray<double>& nodes,
    apf::NewArray<double>& elevatedNodes)
{
  for(int i = 0; i < getNumControlPoints(apf::Mesh::TRIANGLE,P+r); ++i)
    elevatedNodes[i] = 0.;
  raiseBezierTriangle(P,r,nodes,elevatedNodes);
}

template <class T>
static void raiseBezierTet(int P, int r, apf::NewArray<T>& nodes,
    apf::NewArray<T>& elevatedNodes)
{
  for(int i = 0; i <= P+r; ++i){
    for(int j = 0; j <= P+r-i; ++j){
      for(int k = 0; k <= P+r-i-j; ++k){

        for(int l = std::max(0,i-r); l <= std::min(i,P); ++l){
          for(int m = std::max(0,i-l+j-r); m <= std::min(j,P-l); ++m){
            for(int n = std::max(0,i-l+j-m+k-r); n <= std::min(k,P-l-m); ++n){
              elevatedNodes[computeTetNodeIndex(P+r,i,j,k)] +=
                  nodes[computeTetNodeIndex(P,l,m,n)]*quadnomial(P,l,m,n)
                  *quadnomial(r,i-l,j-m,k-n)/quadnomial(P+r,i,j,k);
            }
          }
        }
      }
    }
  }
}

void elevateBezierTet(int P, int r, apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3>& elevatedNodes)
{
  for(int i = 0; i < getNumControlPoints(apf::Mesh::TET,P+r); ++i)
    elevatedNodes[i].zero();
  raiseBezierTet(P,r,nodes,elevatedNodes);
}

static void elevateBezierTetJacobianDet(int P, int r,
    apf::NewArray<double>& nodes,
    apf::NewArray<double>& elevatedNodes)
{
  for(int i = 0; i < getNumControlPoints(apf::Mesh::TET,P+r); ++i)
    elevatedNodes[i] = 0.;
  raiseBezierTet(P,r,nodes,elevatedNodes);
}

typedef void (*ElevateJacobianDetFunction)(int P, int r,
    apf::NewArray<double>& nodes,
    apf::NewArray<double>& elevatedNodes);

const ElevateJacobianDetFunction
elevateBezierJacobianDetArray[apf::Mesh::TYPES] =
{
  NULL,   //vertex
  elevateBezierEdgeJacobianDet,     //edge
  elevateBezierTriangleJacobianDet, //triangle
  NULL,   //quad
  elevateBezierTetJacobianDet,      //tet
  NULL,    //hex
  NULL,    //prism
  NULL     //pyramid
};

void elevateBezierJacobianDet(int type, int P, int r,
    apf::NewArray<double>& nodes,
    apf::NewArray<double>& elevatedNodes)
{
  elevateBezierJacobianDetArray[type](P,r,nodes,elevatedNodes);
}

typedef void (*ElevateFunction)(int P, int r,
    apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3>& elevatedNodes);

const ElevateFunction
elevateBezierArray[apf::Mesh::TYPES] =
{
  NULL,   //vertex
  elevateBezierEdge,     //edge
  elevateBezierTriangle, //triangle
  NULL,   //quad
  elevateBezierTet,      //tet
  NULL,    //hex
  NULL,    //prism
  NULL     //pyramid
};

void elevateBezier(int type, int P, int r,
    apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3>& elevatedNodes)
{
  elevateBezierArray[type](P,r,nodes,elevatedNodes);
}

void elevateBezierCurve(apf::Mesh2* m, apf::MeshEntity* edge, int n, int r)
{
  apf::Element* elem =
      apf::createElement(m->getCoordinateField(),edge);

  apf::Vector3 pt;
  apf::NewArray<apf::Vector3> nodes, elevatedNodes(n+r+1);
  apf::getVectorNodes(elem,nodes);
  elevateBezierEdge(n,r,nodes,elevatedNodes);

  for(int i = 1; i < n+r; ++i)
    m->setPoint(edge,i-1,elevatedNodes[i]);

  apf::destroyElement(elem);
}

}

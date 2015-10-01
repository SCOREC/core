/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crv.h"
#include "crvBezier.h"
#include "crvTables.h"
#include "crvQuality.h"

namespace crv {
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

const ElevateFunction elevateBezierJacobianDet[apf::Mesh::TYPES] =
{
  NULL,   //vertex
  elevateBezierEdgeJacobianDet,     //edge
  elevateBezierTriangleJacobianDet, //triangle
  elevateBezierTetJacobianDet,      //quad
  NULL,      //tet
  NULL,      //hex
  NULL,      //prism
  NULL     //pyramid
};

}

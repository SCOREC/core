/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crv.h"
#include "crvBezier.h"
#include "crvTables.h"

namespace crv {
/*
 * Templating is used for coordinates (Vector3) and det(Jacobian) (double)
 * and is only accessible in this file.
 */
template <class T>
static void raiseBezierEdge(int P, int r, apf::NewArray<T>& nodes,
    apf::NewArray<T>& elevatedNodes)
{
  // re-order nodes, makes life easier
  T temp = nodes[1];
  for (int i = 1; i < P; ++i)
    nodes[i] = nodes[i+1];
  nodes[P] = temp;

  elevatedNodes[0] = nodes[0];
  elevatedNodes[P+r] = nodes[P];
  for(int i = 1; i < P+r; ++i){
    elevatedNodes[i].zero();
    for(int j = std::max(0,i-r); j <= std::min(i,P); ++j)
      elevatedNodes[i] += nodes[j]*binomial(P,j)*binomial(r,i-j)
      /binomial(P+r,i);
  }
}

void elevateBezierEdge(int P, int r, apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3>& elevatedNodes)
{
  raiseBezierEdge(P,r,nodes,elevatedNodes);
}

template <class T>
static void raiseBezierTriangle(int P, int r, apf::NewArray<T>& nodes,
    apf::NewArray<T>& elevatedNodes)
{
  for(int i = 0; i < P+r+1; ++i){
    for(int j = 0; j < P+r+1-i; ++j){
      elevatedNodes[b2[P+r][i][j]].zero();
      for(int k = std::max(0,i-r); k <= std::min(i,P); ++k){
        for(int l = std::max(0,i-k+j-r); l <= std::min(j,P-k); ++l){
          elevatedNodes[b2[P+r][i][j]] += nodes[b2[P][k][l]]
              *trinomial(P,k,l)*trinomial(r,i-k,j-l)/trinomial(P+r,i,j);
        }
      }
    }
  }
}

void elevateBezierTriangle(int P, int r, apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3>& elevatedNodes)
{
  raiseBezierTriangle(P,r,nodes,elevatedNodes);
}

}

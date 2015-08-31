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

/* de Casteljau's algorithm on an edge
 * subNodes[i] corresponds to the i'th edge
 */
template <class T>
static void copyTriangleNodes(int P, apf::NewArray<T>& nodes,
    apf::NewArray<T>& copy)
{
  for (int i = 0; i < (P+1)*(P+2)/2; ++i)
    copy[i] = nodes[i];
}

template <class T>
static void splitEdge(int P, double t, apf::NewArray<T>& nodes,
    apf::NewArray<T> (&subNodes)[2])
{
  // re-order nodes, makes life easier
  T temp = nodes[1];
  for (int i = 1; i < P; ++i)
    nodes[i] = nodes[i+1];
  nodes[P] = temp;

  subNodes[0][0] = nodes[0];
  subNodes[1][P] = nodes[1];
  // go through and find new points,
  // the j'th point on the i'th iteration
  for (int i = 0; i < P; ++i){
    for (int j = 0; j < P-i; ++j){
      nodes[j] = nodes[j]*(1.-t)+nodes[j+1]*t;
    }
    subNodes[0][i+1] = nodes[0];
    subNodes[1][P-i-1] = nodes[P-i-1];
  }
}

void subdivideBezierEdge(int P, double t, apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3> (&subNodes)[2])
{
  splitEdge(P,t,nodes,subNodes);
}

/* de Casteljau's algorithm on a triangle
 * subNodes[i] corresponds to the i'th edge
 */

template <class T, int N>
static void splitTriangle(int P, apf::Vector3& p, apf::NewArray<T>& nodes,
    apf::NewArray<T> (&subNodes)[N], int tri[N])
{
  // set up first two vertices
  for(int t = 0; t < N; ++t)
    subNodes[t][0] = nodes[tri[t]];

  for(int t = 0; t < N; ++t)
    subNodes[t][1] = nodes[(tri[t]+1) % 3];
  // subNodes[t][2] is the split point

  // set up the first edge for all three,
  // using existing edges
  for(int t = 0; t < N; ++t)
    for(int i = 0; i < P-1; ++i)
      subNodes[t][3+i] = nodes[3+tri[t]*(P-1)+i];
  // now each subNodes is filled from 0 to 2+3*(P-1)

  for (int m = 0; m < P; ++m){
    // set up all the nodes for this stage
    for (int i = 0; i < P-m; ++i){
      for (int j = 0; j < P-i-m; ++j){
        int index[3] = {b2[P][i][j],b2[P][i+1][j],b2[P][i][j+1]};

        nodes[index[0]] = nodes[index[0]]*p[0] + nodes[index[1]]*p[1]
                        + nodes[index[2]]*p[2];
      }
    }
    // cycle through the three triangles. each one gets P-m points
    for (int p = 0; p < P-m; ++p){
      int index[3] = {b2[P][P-m-p-1][p],
         b2[P][0][P-m-p-1],b2[P][p][0]};
      for (int t = 0; t < N; ++t)
        subNodes[t][index[0]] = nodes[index[tri[t]]];
    }
  }
}

void subdivideBezierTriangle(int P, apf::Vector3& p,
    apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3> (&subNodes)[3])
{
  int tri[3] = {0,1,2};
  splitTriangle(P,p,nodes,subNodes,tri);
}

/* Four calls of de casteljau's algorithm to subdivide into 4 triangles
 * Uses a non-convex split, which may be unstable, but other work seems
 * to think its okay
 */
void subdivideBezierTriangle(int P, apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3> (&subNodes)[4])
{
  int n = (P+1)*(P+2)/2;
  apf::NewArray<apf::Vector3> tempSubNodes1[1];
  apf::NewArray<apf::Vector3> tempSubNodes2[2];
  tempSubNodes1[0].allocate(n);
  tempSubNodes2[0].allocate(n);
  tempSubNodes2[1].allocate(n);

  int tri1[1] = {1};
  int tri2[2] = {0,1};

  apf::Vector3 p(0.5,0.5,0);
  splitTriangle(P,p,nodes,tempSubNodes2,tri2);
  copyTriangleNodes(P,tempSubNodes2[0],nodes);

  p = apf::Vector3(0,0.5,0.5);
  splitTriangle(P,p,tempSubNodes2[1],tempSubNodes1,tri1);
  copyTriangleNodes(P,tempSubNodes1[0],subNodes[2]);

  tri2[0] = 1; tri2[1] = 2;
  splitTriangle(P,p,nodes,tempSubNodes2,tri2);
  copyTriangleNodes(P,tempSubNodes2[1],subNodes[0]);
  copyTriangleNodes(P,tempSubNodes2[0],nodes);

  p = apf::Vector3(-1,1,1);
  splitTriangle(P,p,nodes,tempSubNodes2,tri2);
  copyTriangleNodes(P,tempSubNodes2[1],subNodes[1]);
  copyTriangleNodes(P,tempSubNodes2[0],subNodes[3]);

}

} // namespace crv

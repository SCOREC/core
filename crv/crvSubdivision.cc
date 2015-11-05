/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crv.h"
#include "crvBezier.h"
#include "crvMath.h"
#include "crvTables.h"
#include "crvQuality.h"

namespace crv {

template <class T>
static void copyTriangleNodes(int P, apf::NewArray<T>& nodes,
    apf::NewArray<T>& copy)
{
  for (int i = 0; i < (P+1)*(P+2)/2; ++i)
    copy[i] = nodes[i];
}

/* de Casteljau's algorithm on an edge
 * subNodes[i] corresponds to the i'th edge
 * P(P+1)/2 additions per split
 */
template <class T>
static void splitEdge(int P, double t, apf::NewArray<T>& nodes,
    apf::NewArray<T> *subNodes)
{
  subNodes[0][0] = nodes[0];
  subNodes[1][P] = nodes[P];
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
  // re-order nodes, makes life easier
  apf::Vector3 temp = nodes[1];
  for (int i = 1; i < P; ++i)
    nodes[i] = nodes[i+1];
  nodes[P] = temp;
  splitEdge(P,t,nodes,subNodes);
}

static void subdivideBezierEdgeJacobianDet(int P, apf::NewArray<double>& nodes,
    apf::NewArray<double> *subNodes)
{
  splitEdge(P,0.5,nodes,subNodes);
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
        unsigned index[3] = {b2[P][i][j],b2[P][i+1][j],b2[P][i][j+1]};

        nodes[index[0]] = nodes[index[0]]*p[0] + nodes[index[1]]*p[1]
                        + nodes[index[2]]*p[2];
      }
    }
    // cycle through the three triangles. each one gets P-m points
    for (int q = 0; q < P-m; ++q){
      unsigned index[3] = {b2[P][P-m-q-1][q],
         b2[P][0][P-m-q-1],b2[P][q][0]};
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
template <class T>
static void splitBezierTriangle(int P, apf::NewArray<T>& nodes,
    apf::NewArray<T> *subNodes)
{
  int n = (P+1)*(P+2)/2;
  apf::NewArray<T> tempSubNodes1[1];
  apf::NewArray<T> tempSubNodes2[2];
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

void subdivideBezierTriangle(int P, apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3> (&subNodes)[4])
{
  splitBezierTriangle(P,nodes,subNodes);
}

static void subdivideBezierTriangleJacobianDet(int P,
    apf::NewArray<double>& nodes,
    apf::NewArray<double> *subNodes)
{
  splitBezierTriangle(P,nodes,subNodes);
}

template <class T, int N>
static void splitTet(int P, apf::Vector3& p, apf::NewArray<T>& nodes,
    apf::NewArray<T> (&subNodes)[N], int tet[N])
{
  // set up first three vertices
  for(int t = 0; t < N; ++t)
  	for (int i = 0; i < 3; ++i)
  		subNodes[t][i] = nodes[apf::tet_tri_verts[tet[t]][i]];

  // subNodes[t][3] is the split point
  // set up edges based on face
  for(int t = 0; t < N; ++t)
  	for(int e = 0; e < 3; ++e)
  		if(!flip_tet_tri_edges[tet[t]][e]) {
  			for(int i = 0; i < P-1; ++i)
  				subNodes[t][4+(P-1)*e+i] = nodes[4+tet_tri_edges[tet[t]][e]*(P-1)+i];
  		} else {
  			for(int i = 0; i < P-1; ++i)
  				subNodes[t][4+(P-1)*e+i] = nodes[4+(tet_tri_edges[tet[t]][e]+1)*(P-1)-1-i];
  		}

  double p3 = 1.-p[0]-p[1]-p[2];
  for (int m = 0; m < P; ++m){
    // set up all the nodes for this stage
    for (int i = 0; i < P-m; ++i){
      for (int j = 0; j < P-i-m; ++j){
      	for (int k = 0; k < P-i-j-m; ++k){
      		unsigned index[4] = {b3[P][i][j][k],b3[P][i+1][j][k],
      				b3[P][i][j+1][k],b3[P][i][j][k+1]};
      		nodes[index[0]] = nodes[index[0]]*p[0] + nodes[index[1]]*p[1]
													+ nodes[index[2]]*p[2] + nodes[index[3]]*p3;
      	}
      }
    }
    // cycle through the tets. each one gets (P-m)*(P-m+1)/2 points
    for (int r = 0; r < P-m; ++r){
    	for (int q = 0; q < P-m-r; ++q){
    		unsigned index[4] = {b3[P][q][r][P-m-1-r-q],b3[P][0][q][r],
    				b3[P][P-m-1-r-q][0][q],b3[P][r][P-m-1-r-q][0]};
    		for (int t = 0; t < N; ++t){
    			subNodes[t][index[0]] = nodes[index[tet[t]]];
    		}
    	}
    }
  }
}

void subdivideBezierTet(int P, apf::Vector3& p,
    apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3> (&subNodes)[4])
{
  int tet[4] = {0,1,2,3};
  splitTet(P,p,nodes,subNodes,tet);
}

void subdivideBezierEntityJacobianDet(int P, int type,
    apf::NewArray<double>& c, apf::NewArray<double>& nodes,
    apf::NewArray<double> *subNodes){
  int typeDim = apf::Mesh::typeDimension[type];
  int n = getNumControlPoints(type,P);
  int numSubdivisions = intpow(2,typeDim);
  for(int k = 0; k < numSubdivisions; ++k){
    for(int j = 0; j < n; ++j){
      subNodes[k][j] = 0.;
      for(int i = 0; i < n; ++i){
        subNodes[k][j] += nodes[i]*c[k*n*n+i+j*n];
      }
    }
  }
}

const SubdivisionFunction subdivideBezierJacobianDet[apf::Mesh::TYPES] =
{
  NULL,   //vertex
  subdivideBezierEdgeJacobianDet,     //edge
  subdivideBezierTriangleJacobianDet, //triangle
  NULL,      //quad
  NULL,      //tet
  NULL,      //hex
  NULL,      //prism
  NULL     //pyramid
};
} // namespace crv

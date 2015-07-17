/******************************************************************************

  Copyright 2015 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/

#ifndef CRVBEZIER_H
#define CRVBEZIER_H

#include "crv.h"
#include "apfDynamicMatrix.h"

namespace crv {

// negative -> flipped relative to canonical
// relies on e0 being always ordered correctly
static int const tet_tri_edges[4][3] =
{{0,1,2},{0,4,3},{1,5,4},{2,5,3}};
static bool const flip_tet_tri_edges[4][3] =
{{0,0,0},{0,0,1},{0,0,1},{1,0,1}};

enum {
  BEZIER,
  GREGORY,
  TYPES
};

// numbers of nodes on
static int const curved_face_internal[2][6] =
{{0,0,1,3,6,10},{0,0,0,6,0,0}};

static int const curved_tet_internal[2][6] =
{{0,0,0,1,4,10},{0,0,0,1,4,10}};

// total numbers of nodes
static int const curved_face_total[2][6] =
{{3,6,10,15,21,28},{0,0,0,18,0,0}};

static int const blended_tet_total[2][6] =
{{4,10,20,34,52,74},{0,0,0,46,0,0}};

static int const curved_tet_total[2][6] =
{{4,10,20,35,56,84},{0,0,0,47,0,0}};

/** \brief shape blending functions */
void BlendedTriangleGetValues(apf::Mesh* m, apf::MeshEntity* e,
    apf::Vector3 const& xi, apf::NewArray<double>& values);

void BlendedTriangleGetLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
    apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& grads);

void BlendedTetGetValues(apf::Mesh* m, apf::MeshEntity* e,
    apf::Vector3 const& xi, apf::NewArray<double>& values);

void BlendedTetGetLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
    apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& grads);

/** \brief get bezier node locations in parameter space */
void getBezierNodeXi(int type, int P, int node, apf::Vector3& xi);

/** \brief compute the matrix to transform between Bezier and Lagrange Points
 *
 \details this is a support function, not actual ever needed.
 Bezier control points, C, can be written as L = A*C
 where A is a matrix of Bernstein polynomials and binomial coefficients.
 To compute the control points from Lagrange points, A^{-1} is used.
 crvBezierPoints.cc contains A^{-1}, precomputed, for
 the nodeXi locations in getBezierNodeXi.

 If, in the future, new nodeXi are defined for interpolating Bezier curves,
 this function can be used to generate the A matrix to invert, as apf
 has no functionality for generic matrix inversion.*/
void getTransformationMatrix(apf::Mesh* m, int type, apf::DynamicMatrix& A);

}

#endif

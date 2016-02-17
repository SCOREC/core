/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRVBEZIER_H
#define CRVBEZIER_H

#include "crv.h"
#include "mth.h"

/** \file crvBezier.h
  * \brief main file for bezier shape functions support */

namespace crv {

/** \brief sets order used in bezier shape functions */
void setOrder(const int order);
/** \brief gets order used in bezier shape functions */
int getOrder();

/** \brief computes node index, use getTriNodeIndex to leverage tables */
int computeTriNodeIndex(int P, int i, int j);
/** \brief computes node index, use getTetNodeIndex to leverage tables */
int computeTetNodeIndex(int P, int i, int j, int k);

/** \brief calculates total number of control points, use
    tables for smaller numbers, this is for quality
    \details This gives the numbers for full bezier shapes,
    this is not accurate for blended shapes */
int getNumControlPoints(int type, int order);
/** \brief calculates number of internal control points, use
    tables for smaller numbers, this is for quality
    \details This gives the numbers for full bezier shapes,
    this is not accurate for blended shapes */
int getNumInternalControlPoints(int type, int order);

/** \brief computes nodes of face f from tet
    \details does not consider alignment, extra care needed */
void getTriNodesFromTetNodes(int f, int P,
    apf::NewArray<apf::Vector3>& tetNodes,
    apf::NewArray<apf::Vector3>& triNodes);
/** \brief computes det(Jacobian) nodes of face f from tet */
void getTriDetJacNodesFromTetDetJacNodes(int f, int P,
    apf::NewArray<double>& tetNodes,
    apf::NewArray<double>& triNodes);
/** \brief gets full set of bezier control points given blended points
    \details this is used for quality assessment of blended shapes,
    and reallocates elemNodes */
void getFullRepFromBlended(apf::Mesh* m, int type,
    apf::NewArray<apf::Vector3>& elemNodes);
/** \brief computes det(Jacobian) for tri from the Bezier conversion
    \details this evaluates an order d(P-1) bezier to compute it
    and is much, much slower than the direct method, but exists for
    comparison */
double computeTriJacobianDetFromBezierFormulation(apf::Mesh* m,
    apf::MeshEntity* e, apf::Vector3& xi);
/** \brief computes det(Jacobian) for tri from the Bezier conversion
    \details this evaluates an order d(P-1) bezier to compute it
    and is much, much slower than the direct method, but exists for
    comparison */
double computeTetJacobianDetFromBezierFormulation(apf::Mesh* m,
    apf::MeshEntity* e, apf::Vector3& xi);

/** \brief shape blending functions
    \details see bezier.tex */
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

/** \brief elevation functions for beziers */
void elevateBezierEdge(int P, int r, apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3>& elevatedNodes);
void elevateBezierTriangle(int P, int r, apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3>& elevatedNodes);
void elevateBezierTet(int P, int r, apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3>& elevatedNodes);
void elevateBezier(int type, int P, int r,
    apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3>& elevatedNodes);

/** \brief subdivision functions for beziers */
void subdivideBezierEdge(int P, double t, apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3> (&subNodes)[2]);
void subdivideBezierTriangle(int P, apf::Vector3& p,
    apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3> (&subNodes)[3]);
void subdivideBezierTriangle(int P, apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3> (&subNodes)[4]);
void subdivideBezierTet(int P, apf::Vector3& p,
    apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3> (&subNodes)[4]);

/** \brief compute the matrix to transform between Bezier and Lagrange Points
    \details this is a support function, not actual ever needed.
    Bezier control points, C, can be written as L = A*C
    where A is a matrix of Bernstein polynomials and binomial coefficients.
    To compute the control points from Lagrange points, A^{-1} is used.
    crvBezierPoints.cc contains A^{-1}, precomputed, for
    the nodeXi locations in getBezierNodeXi. */
void getTransformationMatrix(apf::Mesh* m, apf::MeshEntity* e,
    mth::Matrix<double>& A, const apf::Vector3 *range);

}

#endif

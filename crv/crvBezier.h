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
    \brief main file for support of bezier shape functions */

namespace crv {

/** \brief polynomial part of bernstein polynomial, Bij, Bijk, Bijkl */
double Bij(const int i, const int j,const double u, const double v);
double Bijk(const int i, const int j, const int k, const double u,
    const double v, const double w);
double Bijkl(const int i, const int j, const int k, const int l,
    const double u, const double v, const double w, const double t);

/** \brief a different form of Bij, Bijk, Bijkl */
double Bij(const int ij[], const double xi[]);
double Bijk(const int ijk[], const double xi[]);
double Bijkl(const int ijkl[], const double xi[]);

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

/** \brief Elevate a bezier curve to a higher order
 \details This elevates from nth order to n+rth order
 requires the curve be order n+r in memory already, and
 that the first n points correspond to the lower order curve */
void elevateBezierCurve(apf::Mesh2* m, apf::MeshEntity* edge, int n, int r);

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

/** \brief converts interpolating points to control points
 * \details n is total number of nodes on the shape
 * ne is nodes on the entity, that belong to it
 * c is a coefficient matrix in vector form
 * corresponding to the matrix */
void convertInterpolationPoints(int n, int ne,
    apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<double>& c,
    apf::NewArray<apf::Vector3>& newNodes);

void convertInterpolationPoints(apf::Mesh2* m, apf::MeshEntity* e,
    int n, int ne, apf::NewArray<double>& c);

/** \brief get coefficients for interpolating points to control points
 \details works only for prescribed optimal point locations up to 6th order
 in 2D and */
void getBezierTransformationCoefficients(int P, int type,
    apf::NewArray<double>& c);
void getInternalBezierTransformationCoefficients(apf::Mesh* m, int P, int blend,
    int type, apf::NewArray<double>& c);
void getGregoryTransformationCoefficients(int type, apf::NewArray<double>& c);
void getGregoryBlendedTransformationCoefficients(int blend, int type,
    apf::NewArray<double>& c);

/** \brief a per entity version of above */
void snapToInterpolate(apf::Mesh2* m, apf::MeshEntity* e);

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

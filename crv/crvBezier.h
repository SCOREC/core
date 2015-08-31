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

/** \brief polynomial part of bernstein polynomial */
double Bij(int i, int j,double u, double v);
double Bijk(int i, int j, int k, double u, double v, double w);
double Bijkl(int i, int j, int k, int l,double u,
    double v, double w, double t);

double Bij(int ij[], double xi[]);
double Bijk(int ijk[], double xi[]);
double Bijkl(int ijkl[], double xi[]);

int getTriPointIndex(int P, int i, int j);
int getTetPointIndex(int P, int i, int j, int k);

double computeAlternateTriJacobianDet(apf::Mesh* m,
    apf::MeshEntity* e, apf::Vector3& xi);
double computeAlternateTetJacobianDet(apf::Mesh* m,
    apf::MeshEntity* e, apf::Vector3& xi);
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
/** \brief This is the 2D version. */
void getBezierCurveNodeXi(int type, int P, int node, apf::Vector3& xi);

/** \brief elevation functions for bezier's */
void elevateBezierEdge(int P, int r, apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3>& elevatedNodes);
void elevateBezierTriangle(int P, int r, apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3>& elevatedNodes);

/** \brief subdivision functions for bezier's */
void subdivideBezierEdge(int P, double t, apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3> (&subNodes)[2]);
void subdivideBezierTriangle(int P, apf::Vector3& p,
    apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3> (&subNodes)[3]);
void subdivideBezierTriangle(int P, apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<apf::Vector3> (&subNodes)[4]);

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

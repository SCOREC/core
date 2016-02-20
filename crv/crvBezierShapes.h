/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRVBEZIERSHAPES_H
#define CRVBEZIERSHAPES_H

/** \file crvBezierShapes.h
  * \brief main file for bezier shape functions */

namespace crv {

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

/** \brief typedef for table of shape functions */
typedef void (*bezierShape)(int P,
    apf::Vector3 const& xi,
    apf::NewArray<double>& values);

/** \brief typedef for table of shape function gradients */
typedef void (*bezierShapeGrads)(int P,
    apf::Vector3 const& xi,
    apf::NewArray<apf::Vector3>& grads);

/** \brief table of shape functions */
extern const bezierShape bezier[apf::Mesh::TYPES];
/** \brief table of shape function gradients */
extern const bezierShapeGrads bezierGrads[apf::Mesh::TYPES];

/** \brief Get transformation matrix corresponding to a parametric range
    \details Range is an array of size(num vertices), this is used for
    subdivision, refinement. It is the element transformation matrix,
    see bezier.tex */
void getBezierTransformationMatrix(int type, int P,
    mth::Matrix<double>& A,
    const apf::Vector3 *range);
/** \brief Get transformation matrix of a lower entity in a higher entity
    over a parametric range
    \details Range is an array of size(num vertices), this is used for
    refinement, It is the lower dimensional component, as part of the
    higher array, see bezier.tex */
void getBezierTransformationMatrix(int parentType,
    int childType, int P,
    mth::Matrix<double>& A,
    const apf::Vector3* childRange);

}

#endif

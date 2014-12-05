/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef SPR_H
#define SPR_H

/** \file spr.h
 *  \brief The SPR error estimator interface
 */

#include "apf.h"
#include "apfNew.h"
#include "apfDynamicVector.h"
#include "apfDynamicMatrix.h"

/** \namespace spr
  * \brief All SPR error estimator functions
  */
namespace spr {

/** @brief compute the gradient of a vector or scalar
  *        field at integration points
  * @param f (In) scalar or vector nodal field
  * @param name (In) name of integration point field
  * @param order (In) integration order of accuracy
  */
apf::Field* getGradIPField(apf::Field* f,
                           const char* name, 
                           int order);

/** @brief recover a nodal field using patch recovery
  * @param ip_field (In) integration point field
  */
apf::Field* recoverField(apf::Field* ip_field);

/** @brief run the SPR ZZ error estimator
  * @param f the integration-point input field
  * @param adapt_ratio the fraction of allowable error,
  *                    scales the output size field.
  * @returns a nodal scalar size field with the same distribution
  *          as the mesh coordinate field
  */
apf::Field* getSPRSizeField(apf::Field* f, double adapt_ratio);

/** @brief solve linear least squares problem Ax=b.
  * @param A (In) mxn matrix
  * @param x (Out) nx1 solution vector
  * @param b (In) mx1 right hand side vector
  */
void solveSVD(apf::DynamicMatrix& A,
              apf::DynamicVector& x,
              apf::DynamicVector& b);

/** @brief finds the QR factorization of A
  * @param A the input matrix (rows >= cols)
  * @param V the output representation of Q
  *          as the Householder vectors V(k,:)
  * @param R the output R matrix
  * @returns true iff the matrix A is full rank
  */
bool decompQR(apf::DynamicMatrix& A,
              apf::DynamicMatrix& V,
              apf::DynamicMatrix& R);

/** @brief solves A*x = b given A's QR factorization
  * @param V the Householder vectors from decompQR()
  * @param R the R matrix from decompQR()
  * @param b the right hand side
  * @param x the solution
  */
void solveFromQR(apf::DynamicMatrix& V,
                 apf::DynamicMatrix& R,
                 apf::DynamicVector& b,
                 apf::DynamicVector& x);

}

#endif

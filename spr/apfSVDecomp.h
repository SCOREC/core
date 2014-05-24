/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFSINGULARVALUEDECOMP_H
#define APFSINGULARVALUEDECOMP_H

#include "apfDynamicVector.h"
#include "apfDynamicMatrix.h"

namespace apf {

/** \brief Solves a linear system using singular value decomposition
  * 
  * \details solves the linear system Ax=b using singular value 
  * decomposition, where A can be written:
  * A = U * W * transpose(V) 
  * the solution is found by
  * x = V * inverse(W) * transpose(U) * b
  * since V and W are orthogonal matrices 
  * and inverse(W) is trivial since it is a diagonal matrix 
  *
  * \param A The matrix on which SVD is to be performed and subsequently
  *          inverted
  * \param b The right hand side of the linear system
  * \param maxIters The maximum number of iterations used for SVD
  * \param tol A tolerance for below which singular values are
               assumed to be 0.
  * \param x The vector to be solved for
  */
  void svdSolve(DynamicMatrix& A, DynamicVector& b, 
                int max_iters, double tol, DynamicVector& x);

/** \brief Computes the singular value decomposition of a matrix 
  *
  * \details A = U * W * transpose(V)
  * U - mxn column orthogonal matrix (original matrix A is 
  * updated to U in this routine)
  * w - nxn diagonal matrix whose entries 
  * are called "singular values"
  * (diagonal values W(i,i) are stored in vector w(i))
  * V - nxn orthoganal matrix
  * 
  * \param A The matrix that is decomposed, output as matrix U
  * \param w Vector of singular values
  * \param V Orthoganl matrix
  * \param maxIters Maximum number of iterations for SVD
  */
  void svdDecompose(DynamicMatrix& A, DynamicVector& W,
                    DynamicMatrix& V, int max_iters);
}

#endif

/*
 * Copyright (C) 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_SVD_H
#define APF_SVD_H

#include "apf.h"
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
  void solveSVD(DynamicMatrix& A, 
                DynamicVector& x,
                DynamicVector& b);

}

#endif

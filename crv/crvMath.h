/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRVMATH_H
#define CRVMATH_H

#include "crv.h"
#include <mthQR.h>

/*
 * Brian won't let me put stuff in his mth,
 * because hes a dictator
 * so this is the next best thing
 */

namespace crv {

/** \brief faster power for integers */
inline double intpow(const double b, const int e)
{
  switch (e) {
  case 0: return 1.0;
  case 1: return b;
  case 2: return b*b;
  case 3: return b*b*b;
  case 4: return b*b*b*b;
  case 5: return b*b*b*b*b;
  case 6: return b*b*b*b*b*b;
  default:
    return intpow(b, e-6) * intpow(b, 6);
  }
}

void invertMatrixWithQR(int n, mth::Matrix<double>& A,
    mth::Matrix<double>& Ai);
void invertMatrixWithPLU(int n, mth::Matrix<double>& A,
    mth::Matrix<double>& Ai);
}

#endif


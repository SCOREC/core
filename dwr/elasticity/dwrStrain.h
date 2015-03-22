/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef DWR_STRAIN_H
#define DWR_STRAIN_H

namespace dwr {

template<class Matrix3x3T>
void computeStrain(
    Matrix3x3T const& gradu,
    Matrix3x3T& strain)
{
  zeroMatrix3x3(strain);
  for (int i=0; i < 3; ++i)
  for (int j=0; j < 3; ++j)
    strain[i][j] = 0.5 * ( gradu[i][j] + gradu[j][i] );
}

}

#endif

/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef DWR_STRESS_H
#define DWR_STRESS_H

namespace dwr {

template<class Matrix3x3T>
void computeStress(
    double E,
    double nu,
    Matrix3x3T const& strain,
    Matrix3x3T& stress)
{
  zeroMatrix3x3(stress);
  double lm = ( E * nu ) / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu) );
  double mu = E / ( 2.0 * (1.0 + nu) );
  stress[0][0] = 2.0 * mu * strain[0][0] +
    lm * (strain[0][0] + strain[1][1] + strain[2][2]);
  stress[1][1] = 2.0 * mu * strain[1][1] +
    lm * (strain[0][0] + strain[1][1] + strain[2][2]);
  stress[2][2] = 2.0 * mu * strain[2][2] +
    lm * (strain[0][0] + strain[1][1] + strain[2][2]);
  stress[0][1] = 2.0 * mu * strain[0][1];
  stress[1][2] = 2.0 * mu * strain[1][2];
  stress[2][0] = 2.0 * mu * strain[2][0];
  stress[1][0] = stress[0][1];
  stress[2][1] = stress[1][2];
  stress[0][2] = stress[2][0];
}

}

#endif

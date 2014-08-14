/******************************************************************************

  Copyright 2014 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

  \author Sebastian Rettenberger (rettenbs AT in.tum.de)
*******************************************************************************/

/**\brief Determinant of a 3x3 submatrix of a 4x4 matrix
  * It will always ignore the first column and the row R from 4x4 matrix
  * to get the 3x3 matrix.
  *
  * \tparam R Ignore this row of the 4x4 matrix
  */

#ifndef MA_MATH_H
#define MA_MATH_H

namespace ma {

template<int R>
double determinant3x3(const double m[4][4])
{
	unsigned int x0 = (R > 0 ? 0 : 1);
	unsigned int x1 = (R > 1 ? 1 : 2);
	unsigned int x2 = (R > 2 ? 2 : 3);

	return (m[x0][1]*m[x1][2]*m[x2][3] +
			m[x0][2]*m[x1][3]*m[x2][1] +
			m[x0][3]*m[x1][1]*m[x2][2]) -
		   (m[x2][1]*m[x1][2]*m[x0][3] +
		    m[x2][2]*m[x1][3]*m[x0][1] +
		    m[x2][3]*m[x1][1]*m[x0][2]);
}

/** \brief Determinant of a 4x4 matrix
  */
inline
double determinant4x4(const double m[4][4])
{
	return (m[0][0] * determinant3x3<0>(m) +
			m[2][0] * determinant3x3<2>(m)) -
		   (m[1][0] * determinant3x3<1>(m) +
			m[3][0] * determinant3x3<3>(m));
}

}

#endif // MA_MATH_H

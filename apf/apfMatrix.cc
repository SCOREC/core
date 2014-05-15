/*
 * Copyright (C) 2011-2013 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfMatrix.h"

namespace apf {

Matrix3x3 cross(Vector3 const& u)
{
  return Matrix3x3( 0    ,-u.z(), u.y(),
                    u.z(), 0    ,-u.x(),
                   -u.y(), u.x(), 0    );
}

Matrix3x3 rotate(Vector3 const& u, double a)
{
  Matrix3x3 I(1,0,0,
              0,1,0,
              0,0,1);
  return I*cos(a) + cross(u)*sin(a) + tensorProduct(u,u)*(1-cos(a));
}

/* this is defined outside getFrame because in
   C++98 local types cannot be template arguments,
   and we use the templated std::swap function */
struct SortStruct
{
  int i;
  double m;
};

Matrix3x3 getFrame(Vector3 const& v)
{
  Matrix<3,3> A;
  A[0] = v;
  /* tiny custom code to sort components by absolute value */
  SortStruct s[3] =
  {{0,fabs(v[0])},{1,fabs(v[1])},{2,fabs(v[2])}};
  if (s[2].m > s[1].m)
    std::swap(s[1],s[2]);
  if (s[1].m > s[0].m)
    std::swap(s[0],s[1]);
  if (s[2].m > s[1].m)
    std::swap(s[1],s[2]);
  /* done, components sorted by increasing magnitude */
  int lc = s[0].i;
  int mc = s[1].i;
  int sc = s[2].i;
  /* use the 2D rotation on the largest components
     (rotate v around the */
  A[1][lc] = -v[mc];
  A[1][mc] = v[lc];
  /* and make the last component zero so that A[0] * A[1] == 0 */
  A[1][sc] = 0;
  /* now we have 2 orthogonal (though not unit) vectors, cross
     product gives the third */
  A[2] = cross(A[0],A[1]);
  return A;
}

}

std::ostream& operator<<(std::ostream& s, apf::Matrix3x3 const& v)
{
  s << '\n';
  s << v[0][0] << ' ' << v[0][1] << ' ' << v[0][2] << '\n';
  s << v[1][0] << ' ' << v[1][1] << ' ' << v[1][2] << '\n';
  s << v[2][0] << ' ' << v[2][1] << ' ' << v[2][2] << '\n';
  return s;
}

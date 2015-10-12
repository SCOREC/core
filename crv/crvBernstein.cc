/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crvBezier.h"
#include "crvMath.h"

namespace crv {

double Bij(const int i, const int j, const double u, const double v)
{
  return intpow(u,i)*intpow(v,j);
}

double Bijk(const int i, const int j, const int k, const double u, const double v, const double w)
{
  return intpow(u,i)*intpow(v,j)*intpow(w,k);
}

double Bijkl(const int i, const int j, const int k, const int l,
    const double u, const double v, const double w, const double t)
{
  return intpow(u,i)*intpow(v,j)*intpow(w,k)*intpow(t,l);
}

double Bij(const int ij[], const double xi[])
{
  return Bij(ij[0],ij[1],xi[0],xi[1]);
}

double Bijk(const int ijk[], const double xi[])
{
  return Bijk(ijk[0],ijk[1],ijk[2],xi[0],xi[1],xi[2]);
}

double Bijkl(const int ijkl[], const double xi[])
{
  return Bijkl(ijkl[0],ijkl[1],ijkl[2],ijkl[3],xi[0],xi[1],xi[2],xi[3]);
}

}

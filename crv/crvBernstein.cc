/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include "crvBezier.h"

#include <math.h>

namespace crv {

double Bij(int i, int j, double u, double v)
{
  return pow(u,i)*pow(v,j);
}

double Bijk(int i, int j, int k, double u, double v, double w)
{
  return pow(u,i)*pow(v,j)*pow(w,k);
}

double Bijkl(int i, int j, int k, int l,
    double u, double v, double w, double t)
{
  return pow(u,i)*pow(v,j)*pow(w,k)*pow(t,l);
}

double Bij(int ij[], double xi[])
{
  return Bij(ij[0],ij[1],xi[0],xi[1]);
}

double Bijk(int ijk[], double xi[])
{
  return Bijk(ijk[0],ijk[1],ijk[2],xi[0],xi[1],xi[2]);
}

double Bijkl(int ijkl[], double xi[])
{
  return Bijkl(ijkl[0],ijkl[1],ijkl[2],ijkl[3],xi[0],xi[1],xi[2],xi[3]);
}

}

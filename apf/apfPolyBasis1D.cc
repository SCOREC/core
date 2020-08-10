/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfPolyBasis1D.h"

#include <mth.h>
#include <mth_def.h>
#include <mthQR.h>
#include <pcu_util.h>
#include <PCU.h>

namespace apf {

static void getPoints(
    int order,
    apf::NewArray<double>& points,
    int type)
{
  int np = order + 1;
  points.allocated() ? points.resize(np) : points.allocate(np);
  switch (type)
  {
    case GAUSS_LEGENDRE:
    {
      getGaussLegendrePoints(np, &points[0]);
      break;
    }
    case GAUSS_LOBATTO:
    {
      getGaussLobattoPoints(np, &points[0]);
      break;
    }
    default:
      PCU_ALWAYS_ASSERT_VERBOSE(false,
      	  "type should be either GAUSS_LEGENDRE or GAUSS_LOBATTO!");
  }
}

void getGaussLegendrePoints(int np, double* pts)
{
  switch (np)
  {
    case 1:
      pts[0] = 0.5;
      return;
    case 2:
      pts[0] = 0.21132486540518711775;
      pts[1] = 0.78867513459481288225;
      return;
    case 3:
      pts[0] = 0.11270166537925831148;
      pts[1] = 0.5;
      pts[2] = 0.88729833462074168852;
      return;
  }

  const int n = np;
  const int m = (n+1)/2;

  for (int i = 1; i <= m; i++)
  {
     double z = cos(M_PI * (i - 0.25) / (n + 0.5));
     double pp, p1, dz, xi = 0.;
     bool done = false;
     while (1)
     {
        double p2 = 1;
        p1 = z;
        for (int j = 2; j <= n; j++)
        {
           double p3 = p2;
           p2 = p1;
           p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j;
        }
        // p1 is Legendre polynomial
        pp = n * (z*p1-p2) / (z*z - 1);
        if (done) { break; }
        dz = p1/pp;
        if (fabs(dz) < 1e-16)
        {
           done = true;
           // map the new point (z-dz) to (0,1):
           xi = ((1 - z) + dz)/2; // (1 - (z - dz))/2 has bad round-off
           // continue the computation: get pp at the new point, then exit
        }
        // update: z = z - dz
        z -= dz;
     }
     pts[i-1] = xi;
     pts[n-i] = 1 - xi;
  }
}

void getGaussLobattoPoints(int /*np*/, double* /*pts*/)
{ /* implement Gauss Lobatto points. Later when needed. */
}


void getOpenPoints(
    int order,
    apf::NewArray<double>& op,
    int type)
{
  getPoints(order, op, type);
}

void getClosedPoints(
    int order,
    apf::NewArray<double>& cp,
    int type)
{
    getPoints(order, cp, type);
}

void getChebyshevT(int order, double xi, double* u)
{
  // recursive definition, z in [-1,1]
  // T_0(z) = 1,  T_1(z) = z
  // T_{n+1}(z) = 2*z*T_n(z) - T_{n-1}(z)
  double z;
  u[0] = 1.;
  if (order == 0) { return; }
  u[1] = z = 2.*xi - 1.;
  for (int n = 1; n < order; n++)
  {
     u[n+1] = 2*z*u[n] - u[n-1];
  }
}

void getChebyshevT(int order, double xi, double* u, double* d)
{
  // recursive definition, z in [-1,1]
  // T_0(z) = 1,  T_1(z) = z
  // T_{n+1}(z) = 2*z*T_n(z) - T_{n-1}(z)
  // T'_n(z) = n*U_{n-1}(z)
  // U_0(z) = 1  U_1(z) = 2*z
  // U_{n+1}(z) = 2*z*U_n(z) - U_{n-1}(z)
  // U_n(z) = z*U_{n-1}(z) + T_n(z) = z*T'_n(z)/n + T_n(z)
  // T'_{n+1}(z) = (n + 1)*(z*T'_n(z)/n + T_n(z))
  double z;
  u[0] = 1.;
  d[0] = 0.;
  if (order == 0) { return; }
  u[1] = z = 2.*xi - 1.;
  d[1] = 2.;
  for (int n = 1; n < order; n++)
  {
     u[n+1] = 2*z*u[n] - u[n-1];
     d[n+1] = (n + 1)*(z*d[n]/n + 2*u[n]);
  }
}

void getChebyshevT(int order, double xi, double* u, double* d, double* dd)
{
  // recursive definition, z in [-1,1]
  // T_0(z) = 1,  T_1(z) = z
  // T_{n+1}(z) = 2*z*T_n(z) - T_{n-1}(z)
  // T'_n(z) = n*U_{n-1}(z)
  // U_0(z) = 1  U_1(z) = 2*z
  // U_{n+1}(z) = 2*z*U_n(z) - U_{n-1}(z)
  // U_n(z) = z*U_{n-1}(z) + T_n(z) = z*T'_n(z)/n + T_n(z)
  // T'_{n+1}(z) = (n + 1)*(z*T'_n(z)/n + T_n(z))
  // T''_{n+1}(z) = (n + 1)*(2*(n + 1)*T'_n(z) + z*T''_n(z)) / n
  double z;
  u[0] = 1.;
  d[0] = 0.;
  dd[0]= 0.;
  if (order == 0) { return; }
  u[1] = z = 2.*xi - 1.;
  d[1] = 2.;
  dd[1] = 0;
  for (int n = 1; n < order; n++)
  {
     u[n+1] = 2*z*u[n] - u[n-1];
     d[n+1] = (n + 1)*(z*d[n]/n + 2*u[n]);
     dd[n+1] = (n + 1)*(2.*(n + 1)*d[n] + z*dd[n])/n;
  }
}

}

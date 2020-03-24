/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfShape.h"
#include "apfMesh.h"
#include "apfFieldOf.h"
#include "apfElement.h"
#include "apfVectorElement.h"
#include "mth.h"
#include <pcu_util.h>
#include <mthQR.h>
#include <PCU.h>

#include <iostream>

namespace apf {

enum {
  GAUSS_LEGENDRE,
  GAUSS_LOBATTO
};

static void getGaussLegendrePoints(int np, double* pts)
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

static void getGaussLobattoPoints(int /*np*/, double* /*pts*/)
{ /* TODO implement Gauss Lobatto points. Later when needed. */
};

static const double* getPoints(int order, const int type)
{
  int np = order + 1;
  double* points = new double[np];
  switch (type)
  {
    case GAUSS_LEGENDRE:
    {
      getGaussLegendrePoints(np, points);
      break;
    }
    case GAUSS_LOBATTO:
    {
        getGaussLobattoPoints(np, points);
      break;
    }
  }
  return points;
}

static const double* getOpenPoints(int order, const int type = GAUSS_LEGENDRE)
{
  return getPoints(order, type);
}

static const double* getClosedPoints(int order, const int type = GAUSS_LOBATTO)
{
    return getPoints(order, type);
}

static void getChebyshevT(int order, double xi, double* u)
{
  // TODO implement Chebyshev
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

static void getChebyshevT(int order, double xi, double* u, double* d)
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

static void getChebyshevT(int order, double xi, double* u, double* d, double* dd)
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

class Nedelec: public FieldShape {
  public:
    const char* getName() const { return "Nedelec"; }
    Nedelec(int order) : P(order) {
    }
    class Vertex : public apf::EntityShape
    {
    public:
      void getValues(apf::Mesh*, apf::MeshEntity*,
	  apf::Vector3 const&, apf::NewArray<double>& values) const
      {
	(void)values;
	// TODO inform the user that this is not implemented and abort()
      }
      void getLocalGradients(apf::Mesh*, apf::MeshEntity*,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      }
      int countNodes() const {return 0;}
      void alignSharedNodes(apf::Mesh*,
	  apf::MeshEntity*, apf::MeshEntity*, int order[])
      {
	(void)order;
      }
    };
    class Edge : public apf::EntityShape
    {
    public:
      void getValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<double>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getValues not implemented for Nedelec Edges. Aborting()!");
      }
      void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getLocalGradients not implemented for Nedelec Edges. Aborting()!");
      }
      int countNodes() const {return 0; /* TODO update this*/}
      void alignSharedNodes(apf::Mesh*,
	  apf::MeshEntity*, apf::MeshEntity*, int order[])
      {
	(void)order;
      }
      void getVectorValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	// TODO: to be completed
      }
      void getLocalVectorCurls(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	// TODO: to be completed
      }
    };
    class Triangle : public apf::EntityShape
    {
    public:
      int getOrder() {return P;}
      void getValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<double>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getValues not implemented for Nedelec Triangle. Try getVectorValues. Aborting()!");
      }
      void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getLocalGradients not implemented for Nedelec Triangle. Aborting()!");
      }
      int countNodes() const {return P*(P+2);}
      void alignSharedNodes(apf::Mesh*,
	  apf::MeshEntity*, apf::MeshEntity*, int order[])
      {
	(void)order;
      }
      void getVectorValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& shapes) const
      {
      	const int pm1 = P - 1;
        const int p = P;

      	apf::NewArray<double> shape_x(p);
      	apf::NewArray<double> shape_y(p);
      	apf::NewArray<double> shape_l(p);

        int dof = countNodes();
        mth::Matrix<double> u(dof, dim);

        double x = xi[0]; double y = xi[1];

        getChebyshevT(pm1, x, &shape_x[0]);
        getChebyshevT(pm1, y, &shape_y[0]);
        getChebyshevT(pm1, 1. - x - y, &shape_l[0]);

        int n = 0;
        for (int j = 0; j <= pm1; j++)
          for (int i = 0; i + j <= pm1; i++)
          {
            double s = shape_x[i]*shape_y[j]*shape_l[pm1-i-j];
            u(n,0) = s;  u(n,1) = 0.;  n++;
            u(n,0) = 0.;  u(n,1) = s;  n++;
          }
        for (int j = 0; j <= pm1; j++)
        {
          double s = shape_x[pm1-j]*shape_y[j];
          u(n,0) = s*(y - c);  u(n,1) = -s*(x - c); n++;
        }


      	computeTi();
        mth::Matrix<double> S(dof, dim);
      	for(int i = 0; i < dim; i++) // S = Ti * u
        {
          apf::Vector<dof> B;
          apf::Vector<double> X(dof);
          for (int j = 0; j < dof; j++) B[j] = u(i,j); // populate b
          mth::solveFromQR(Q, R, B, X);
          for (int j = 0; j < dof; j++)  S(i,j) = X[j]; // populate S with x
        }

        shapes.allocate(dof);
        for (int i = 0; i < dof; i++) // populate y
        {
      	  shapes(i) = apf::Vector3( S(i,0), S(i,1), 0.0 );
        }
      }
      void getLocalVectorCurls(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const& xi, apf::NewArray<apf::Vector3>&) const
      {
      	const int pm1 = P - 1;
        const int p = P;

      	apf::NewArray<double> shape_x(p);
      	apf::NewArray<double> shape_y(p);
      	apf::NewArray<double> shape_l(p);
      	apf::NewArray<double> dshape_x(p);
      	apf::NewArray<double> dshape_y(p);
      	apf::NewArray<double> dshape_l(p);

      	int dof = countNodes();
        apf::Vector<double> curlu(dof);

        double x = xi[0]; double y = xi[1];

        getChebyshevT(pm1, x, &shape_x[0], &dshape_x[0]);
        getChebyshevT(pm1, y, &shape_y[0], &dshape_y[0]);
        getChebyshevT(pm1, 1. - x - y, &shape_l[0], &dshape_l[0]);

        int n = 0;
        for (int j = 0; j <= pm1; j++)
          for (int i = 0; i + j <= pm1; i++)
          {
            int l = pm1-i-j;
            const double dx = (dshape_x[i]*shape_l[l] -
                          shape_x(i)*dshape_l(l)) * shape_y(j);
            const double dy = (dshape_y[j]*shape_l[l] -
                          shape_y[j]*dshape_l[l]) * shape_x[i];

            curlu(n++) = -dy;
            curlu(n++) =  dx;
          }
        for (int j = 0; j <= pm1; j++)
        {
          int i = pm1 - j;
          // curl of shape_x(i)*shape_y(j) * (ip.y - c, -(ip.x - c), 0):
          curlu(n++) = -((dshape_x[i]*(x - c) + shape_x[i]) * shape_y[j] +
                     (dshape_y[j]*(y - c) + shape_y[j]) * shape_x[i]);
        }

        computeTi();
        apf::Vector<double> X(dof);
        mth::solveFromQR(Q, R, curlu, X);
      }
    private:
      int P; // polyonmial order
      int dim = 2; // reference element dim
      double c = 1./3.; // center of tri
      const double tk[8] = {1.,0.,  -1.,1.,  0.,-1.,  0.,1.};
      mth::Matrix<double> Q;
      mth::Matrix<double> R;
      void computeTi()
      {
     	  int non = countNodes();

        const double *eop = getOpenPoints(P - 1);
        const double *iop = (P > 1) ? getOpenPoints(P - 2) : NULL;

        const int p = P, pm1 = P - 1, pm2 = P - 2;

      	apf::NewArray<double> shape_x(p);
      	apf::NewArray<double> shape_y(p);
      	apf::NewArray<double> shape_l(p);

       	apf::DynamicArray<apf::Vector3> nodes (non);
      	apf::DynamicArray<int> dof2tk (non);

        int o = 0;
      	// Edge loops to get nodes and dof2tk for edges
        for (int i = 0; i < P; i++)  // (0,1)
        {
          nodes[o][0] = eop[i];  nodes[o][1] = 0.;  nodes[o][2] = 0.;
          dof2tk[o++] = 0;
        }
        for (int i = 0; i < P; i++)  // (1,2)
        {
          nodes[o][0] = eop[pm1-i];  nodes[o][1] = eop[i];  nodes[o][2] = 0.;
          dof2tk[o++] = 1;
        }
        for (int i = 0; i < P; i++)  // (2,0)
        {
          nodes[o][0] = 0.;  nodes[o][1] = eop[pm1-i];  nodes[o][2] = 0.;
          dof2tk[o++] = 2;
        }

      	// Region loops to get nodes and dof2tk for regions
        for (int j = 0; j <= pm2; j++)
          for (int i = 0; i + j <= pm2; i++)
              {
                double w = iop[i] + iop[j] + iop[pm2-i-j];
                nodes[o][0] = iop[i]/w;  nodes[o][1] = iop[j]/w;  nodes[o][2] = 0.;
                dof2tk[o++] = 0;
                nodes[o][0] = iop[i]/w;  nodes[o][1] = iop[j]/w;  nodes[o][2] = 0.;
                dof2tk[o++] = 3;
              }

      	// Populate T
      	mth::Matrix<double> T(non,non); // T(i,j)
        for (int m = 0; m < non; m++)
        {
          const double *tm = tk + 2*dof2tk[m];
          o = 0;

          double x = nodes[o][0]; double y = nodes[o][1];

          getChebyshevT(pm1, x, &shape_x[0]);
          getChebyshevT(pm1, y, &shape_y[0]);
          getChebyshevT(pm1, 1. - x - y, &shape_l[0]);

          for (int j = 0; j <= pm1; j++)
            for (int i = 0; i + j <= pm1; i++)
            {
              double s = shape_x[i]*shape_y[j]*shape_l[pm1-i-j];
              T(o++, m) = s * tm[0];
              T(o++, m) = s * tm[1];
            }
          for (int j = 0; j <= pm1; j++)
          {
            T(o++, m) =
              shape_x[pm1-j]*shape_y[j]*((y - c)*tm[0] - (x - c)*tm[1]);
          }
        }

        unsigned rank = mth::decomposeQR(T, Q, R);
      }
    };
    class Tetrahedron : public apf::EntityShape
    {
    public:
      int getOrder() {return P;}
      void getValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<double>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getValues not implemented for Nedelec Tetrahedron. Try getVectorValues. Aborting()!");
      }
      void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getLocalGradients not implemented for Nedelec Tetrahedron. Aborting()!");
      }
      int countNodes() const {return (P+3)*(P+2)*P/2;}
      void alignSharedNodes(apf::Mesh*,
	  apf::MeshEntity*, apf::MeshEntity*, int order[])
      {
	(void)order;
      }
      void getVectorValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& shapes) const
      {
        const int pm1 = P - 1;
        const int p = P;

      	apf::NewArray<double> shape_x(p);
      	apf::NewArray<double> shape_y(p);
      	apf::NewArray<double> shape_z(p);
      	apf::NewArray<double> shape_l(p);

        int dof = countNodes();
        mth::Matrix<double> u(dof, dim);

        double x = xi[0]; double y = xi[1]; double z = xi[2];

        getChebyshevT(pm1, x, &shape_x[0]);
        getChebyshevT(pm1, y, &shape_y[0]);
        getChebyshevT(pm1, z, &shape_z[0]);
        getChebyshevT(pm1, 1. - x - y - z, &shape_l[0]);

        int n = 0;
        for (int k = 0; k <= pm1; k++)
          for (int j = 0; j + k <= pm1; j++)
            for (int i = 0; i + j + k <= pm1; i++)
            {
              double s = shape_x[i]*shape_y[j]*shape_z[k]*shape_l[pm1-i-j-k];
                    u(n,0) =  s;  u(n,1) = 0.;  u(n,2) = 0.;  n++;
                    u(n,0) = 0.;  u(n,1) =  s;  u(n,2) = 0.;  n++;
                    u(n,0) = 0.;  u(n,1) = 0.;  u(n,2) =  s;  n++;
                 }
        for (int k = 0; k <= pm1; k++)
          for (int j = 0; j + k <= pm1; j++)
          {
            double s = shape_x[pm1-j-k]*shape_y[j]*shape_z[k];
            u(n,0) = s*(y - c);  u(n,1) = -s*(x - c);  u(n,2) =  0.;  n++;
            u(n,0) = s*(z - c);  u(n,1) =  0.;  u(n,2) = -s*(x - c);  n++;
          }
        for (int k = 0; k <= pm1; k++)
        {
          double s = shape_y[pm1-k]*shape_z[k];
          u(n,0) = 0.;  u(n,1) = s*(z - c);  u(n,2) = -s*(y - c);  n++;
        }

      	computeTi();
        mth::Matrix<double> S(dof, dim);
      	for(int i = 0; i < dim; i++) // S = Ti * u
        {
          apf::Vector<double> B(dof);
          apf::Vector<double> X(dof);
          for (int j = 0; j < dof; j++) B[j] = u(i,j); // populate b
          mth::solveFromQR(Q, R, B, X);
          for (int j = 0; j < dof; j++)  S(i,j) = X[j]; // populate S with x
        }

        shapes.allocate(dof);
        for (int i = 0; i < dof; i++) // populate y
        {
      	  shapes(i) = apf::Vector3( S(i,0), S(i,1), S(i,2) );
        }
      }
      void getLocalVectorCurls(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& curl_shapes) const
      {
      	const int pm1 = P - 1;
        const int p = P;

      	apf::NewArray<double> shape_x(p);
      	apf::NewArray<double> shape_y(p);
      	apf::NewArray<double> shape_z(p);
      	apf::NewArray<double> shape_l(p);
      	apf::NewArray<double> dshape_x(p);
      	apf::NewArray<double> dshape_y(p);
      	apf::NewArray<double> dshape_z(p);
      	apf::NewArray<double> dshape_l(p);

      	int dof = countNodes();
        mth::Matrix<double> u(dof, dim);

        double x = xi[0]; double y = xi[1]; double z = xi[2];

        getChebyshevT(pm1, x, &shape_x[0], &dshape_x[0]);
        getChebyshevT(pm1, y, &shape_y[0], &dshape_y[0]);
        getChebyshevT(pm1, z, &shape_z[0], &dshape_z[0]);
        getChebyshevT(pm1, 1. - x - y - z, &shape_l[0], &dshape_l[0]);

        int n = 0;
        for (int k = 0; k <= pm1; k++)
          for (int j = 0; j + k <= pm1; j++)
            for (int i = 0; i + j + k <= pm1; i++)
            {
              int l = pm1-i-j-k;
              const double dx = (dshape_x[i]*shape_l[l] -
                             shape_x[i]*dshape_l[l])*shape_y[j]*shape_z[k];
              const double dy = (dshape_y[j]*shape_l[l] -
                             shape_y[j]*dshape_l[l])*shape_x[i]*shape_z[k];
              const double dz = (dshape_z[k]*shape_l[l] -
                             shape_z[k]*dshape_l[l])*shape_x[i]*shape_y[j];

              u(n,0) =  0.;  u(n,1) =  dz;  u(n,2) = -dy;  n++;
              u(n,0) = -dz;  u(n,1) =  0.;  u(n,2) =  dx;  n++;
              u(n,0) =  dy;  u(n,1) = -dx;  u(n,2) =  0.;  n++;
            }
        for (int k = 0; k <= pm1; k++)
          for (int j = 0; j + k <= pm1; j++)
          {
            int i = pm1 - j - k;
            // s = shape_x(i)*shape_y(j)*shape_z(k);
            // curl of s*(ip.y - c, -(ip.x - c), 0):
            u(n,0) =  shape_x[i]*(x - c)*shape_y[j]*dshape_z[k];
            u(n,1) =  shape_x[i]*shape_y[j]*(y - c)*dshape_z[k];
            u(n,2) =  -((dshape_x[i]*(x - c) + shape_x[i])*shape_y[j]*shape_z[k] +
                      (dshape_y(j)*(y - c) + shape_y[j])*shape_x[i]*shape_z[k]);
            n++;
            // curl of s*(ip.z - c, 0, -(ip.x - c)):
            u(n,0) = -shape_x[i]*(x - c)*dshape_y[j]*shape_z[k];
            u(n,1) = (shape_x[i]*shape_y[j]*(dshape_z[k]*(z - c) + shape_z[k]) +
                     (dshape_x[i]*(x - c) + shape_x[i])*shape_y[j]*shape_z[k]);
            u(n,2) = -shape_x[i]*dshape_y[j]*shape_z[k]*(z - c);
            n++;
          }
        for (int k = 0; k <= pm1; k++)
        {
          int j = pm1 - k;
          // curl of shape_y(j)*shape_z(k)*(0, ip.z - c, -(ip.y - c)):
          u(n,0) = -((dshape_y[j]*(y - c) + shape_y[j])*shape_z[k] +
                   shape_y[j]*(dshape_z[k]*(z - c) + shape_z[k]));
          u(n,1) = 0.;
          u(n,2) = 0.;  n++;
        }

        computeTi();
        mth::Matrix<double> S(dof, dim);
      	for(int i = 0; i < dim; i++) // S = Ti * u
        {
          apf::Vector<dof> B;
          apf::Vector<double> X(dof);
          for (int j = 0; j < dof; j++) B[j] = u(i,j); // populate b
          mth::solveFromQR(Q, R, B, X);
          for (int j = 0; j < dof; j++)  S(i,j) = X[j]; // populate S with x
        }

        curl_shapes.allocate(dof);
        for (int i = 0; i < dof; i++) // populate y
        {
      	  curl_shapes(i) = apf::Vector3( S(i,0), S(i,1), S(i,2) );
        }
      }
    private:
      int P;
      int dim = 3;
      double c = 1./4.;
      const double tk[18] = {1.,0.,0.,  -1.,1.,0.,  0.,-1.,0.,  0.,0.,1.,  -1.,0.,1.,  0.,-1.,1.}; // edge directions
      mth::Matrix<double> Q;
      mth::Matrix<double> R;
      void computeTi()
      {
     	  int non = countNodes();

        const double *eop = getOpenPoints(P - 1);
        const double *fop = (P > 1) ? getOpenPoints(P - 2) : NULL;
        const double *iop = (P > 2) ? getOpenPoints(P - 3) : NULL;

        const int p = P, pm1 = P - 1, pm2 = P - 2, pm3 = P - 3;

      	apf::NewArray<double> shape_x(p);
      	apf::NewArray<double> shape_y(p);
      	apf::NewArray<double> shape_z(p);
      	apf::NewArray<double> shape_l(p);

       	apf::DynamicArray<apf::Vector3> nodes (non);
      	apf::DynamicArray<int> dof2tk (non);
	      nodes.setSize(non);
      	dof2tk.setSize(non);

        int o = 0;
      	// Edge loops to get nodes and dof2tk for edges
        for (int i = 0; i < P; i++)  // (0,1)
        {
          nodes[o][0] = eop[i];  nodes[o][1] = 0.;  nodes[o][2] = 0.;
          dof2tk[o++] = 0;
        }
        for (int i = 0; i < P; i++)  // (1,2)
        {
          nodes[o][0] = eop[pm1-i];  nodes[o][1] = eop[i];  nodes[o][2] = 0.;
          dof2tk[o++] = 1;
        }
        for (int i = 0; i < P; i++)  // (2,0)
        {
          nodes[o][0] = 0.;  nodes[o][1] = eop[pm1-i];  nodes[o][2] = 0.;
          dof2tk[o++] = 2;
        }
        for (int i = 0; i < P; i++)  // (0,3)
        {
          nodes[o][0] = 0.;  nodes[o][1] = 0.;  nodes[o][2] = eop[i];
          dof2tk[o++] = 3;
        }
        for (int i = 0; i < P; i++)  // (1,3)
        {
          nodes[o][0] = eop[pm1-i];  nodes[o][1] = 0.;  nodes[o][2] = eop[i];
          dof2tk[o++] = 4;
        }
        for (int i = 0; i < P; i++)  // (2,3)
        {
          nodes[o][0] = 0.;  nodes[o][1] = eop[pm1-i];  nodes[o][2] = eop[i];
          dof2tk[o++] = 5;
        }

      	// Face loops to get nodes and dof2tk for faces
        for (int j = 0; j <= pm2; j++)  // (0,1,2)
          for (int i = 0; i + j <= pm2; i++)
          {
            double w = fop[i] + fop[j] + fop[pm2-i-j];
            nodes[o][0] = fop[i]/w;  nodes[o][1] = fop[j]/w;  nodes[o][2] = 0.;
            dof2tk[o++] = 0;
            nodes[o][0] = fop[i]/w;  nodes[o][1] = fop[j]/w;  nodes[o][2] = 0.;
            dof2tk[o++] = 2;
          }
        for (int j = 0; j <= pm2; j++)  // (0,1,3)
          for (int i = 0; i + j <= pm2; i++)
          {
            double w = fop[i] + fop[j] + fop[pm2-i-j];
            nodes[o][0] = fop[i]/w;  nodes[o][1] = 0.;  nodes[o][2] = fop[j]/w;
            dof2tk[o++] = 0;
            nodes[o][0] = fop[i]/w;  nodes[o][1] = 0.;  nodes[o][2] = fop[j]/w;
            dof2tk[o++] = 3;
          }
        for (int j = 0; j <= pm2; j++)  // (1,2,3)
          for (int i = 0; i + j <= pm2; i++)
          {
            double w = fop[i] + fop[j] + fop[pm2-i-j];
            nodes[o][0] = fop[pm2-i-j]/w;  nodes[o][1] = fop[i]/w;  nodes[o][2] = fop[j]/w;
            dof2tk[o++] = 1;
            nodes[o][0] = fop[pm2-i-j]/w;  nodes[o][1] = fop[i]/w;  nodes[o][2] = fop[j]/w;
            dof2tk[o++] = 4;
          }
        for (int j = 0; j <= pm2; j++)  // (0,2,3)
          for (int i = 0; i + j <= pm2; i++)
          {
            double w = fop[i] + fop[j] + fop[pm2-i-j];
            nodes[o][0] = 0.;  nodes[o][1] = fop[i]/w;  nodes[o][2] = fop[j]/w;
            dof2tk[o++] = 2;
            nodes[o][0] = 0.;  nodes[o][1] = fop[i]/w;  nodes[o][2] = fop[j]/w;
            dof2tk[o++] = 3;
          }

      	// Region loops to get nodes and dof2tk for regions
        for (int k = 0; k <= pm3; k++)
          for (int j = 0; j + k <= pm3; j++)
            for (int i = 0; i + j + k <= pm3; i++)
              {
                double w = iop[i] + iop[j] + iop[k] + iop[pm3-i-j-k];
                nodes[o][0] = iop[i]/w;  nodes[o][1] = iop[j]/w;  nodes[o][2] = iop[k]/w;
                dof2tk[o++] = 0;
                nodes[o][0] = iop[i]/w;  nodes[o][1] = iop[j]/w;  nodes[o][2] = iop[k]/w;
                dof2tk[o++] = 2;
                nodes[o][0] = iop[i]/w;  nodes[o][1] = iop[j]/w;  nodes[o][2] = iop[k]/w;
                dof2tk[o++] = 3;
              }

      	// Populate T
      	mth::Matrix<double> T(non,non); // T(i,j)
        for (int m = 0; m < non; m++)
        {
          const double *tm = tk + 3*dof2tk[m];
          o = 0;

          double x = nodes[o][0]; double y = nodes[o][1]; double z = nodes[o][2];

          getChebyshevT(pm1, x, &shape_x[0]);
          getChebyshevT(pm1, y, &shape_y[0]);
          getChebyshevT(pm1, z, &shape_z[0]);
          getChebyshevT(pm1, 1. - x - y - z, &shape_l[0]);

          for (int k = 0; k <= pm1; k++)
            for (int j = 0; j + k <= pm1; j++)
              for (int i = 0; i + j + k <= pm1; i++)
              {
                double s = shape_x[i]*shape_y[j]*shape_z[k]*shape_l[pm1-i-j-k];
                T(o++, m) = s * tm[0];
                T(o++, m) = s * tm[1];
                T(o++, m) = s * tm[2];
              }
          for (int k = 0; k <= pm1; k++)
            for (int j = 0; j + k <= pm1; j++)
              {
                double s = shape_x[pm1-j-k]*shape_y[j]*shape_z[k];
                T(o++, m) = s*((y - c)*tm[0] - (x - c)*tm[1]);
                T(o++, m) = s*((z - c)*tm[0] - (x - c)*tm[2]);
              }
          for (int k = 0; k <= pm1; k++)
          {
            T(o++, m) =
              shape_y[pm1-k]*shape_z[k]*((z - c)*tm[1] - (y - c)*tm[2]);
          }
        }

        unsigned rank = mth::decomposeQR(T, Q, R);
      }
    };
};

};

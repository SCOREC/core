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
#include "apfPolyBasis1D.h"
#include <mth.h>
#include <mth_def.h>
#include <mthQR.h>
#include <pcu_util.h>
#include <PCU.h>

#include <iostream>
using namespace std;

namespace apf {

static void alignFaceNodes(
    int init_order[],
    int final_order[],
    int nodes, int r, bool f)
{
  if (r == 0 && !f)     // CASE 1
  {
    int ind = 0;
    int offset = 0;
    for (int i = 0; i < nodes; i+=2)
    {
      final_order[i] = init_order[offset];
      final_order[i+1] = init_order[offset];
      offset += 2;
    }
    ind += nodes;
    offset = 1;
    for (int i = 0 ; i < nodes; i+=2)
    {
      final_order[i+ind] = init_order[offset];
      final_order[i+1+ind] = init_order[offset];
      offset += 2;
    }
    ind += nodes;
    for (int i = 0; i < nodes; i++)
      final_order[i+ind] = (i%2) ? 0 : 1;
    ind += nodes;
    for (int i = 0; i < nodes; i++)
      final_order[i+ind] = (i%2) ? 1 : 0;
  }
  else if (r == 1 && !f) // CASE 2
  {
    int ind = 0;
    int offset = 0;
    for (int i = 0; i < nodes; i+=2)
    {
      final_order[i] = -init_order[offset]-1;
      final_order[i+1] = -init_order[offset]-1;
      offset += 2;
    }
    ind += nodes;
    offset = 1;
    for (int i = 0 ; i < nodes; i+=2)
    {
      final_order[i+ind] = init_order[offset];
      final_order[i+1+ind] = init_order[offset];
      offset += 2;
    }
    ind += nodes;
    for (int i = 0; i < nodes; i++)
      final_order[i+ind] = 1;
    ind += nodes;
    for (int i = 0; i < nodes; i++)
      final_order[i+ind] = (i%2) ? 0 : 1;
  }
  else if (r == 2 && !f) // CASE 3
  {
    int ind = 0;
    int offset = 0;
    for (int i = 0; i < nodes; i+=2)
    {
      final_order[i] = init_order[offset];
      final_order[i+1] = init_order[offset];
      offset += 2;
    }
    ind += nodes;
    offset = 1;
    for (int i = 0 ; i < nodes; i+=2)
    {
      final_order[i+ind] = -init_order[offset]-1;
      final_order[i+1+ind] = -init_order[offset]-1;
      offset += 2;
    }
    ind += nodes;
    for (int i = 0; i < nodes; i++)
      final_order[i+ind] = (i%2) ? 1 : 0;
    ind += nodes;
    for (int i = 0; i < nodes; i++)
      final_order[i+ind] = 1;
  }
  else if (r == 0 && f)  // CASE 4
  {
    int ind = 0;
    int offset = 0;
    for (int i = 0; i < nodes; i+=2)
    {
      final_order[i] = init_order[offset];
      final_order[i+1] = init_order[offset];
      offset += 2;
    }
    ind += nodes;
    offset = 1;
    for (int i = 0 ; i < nodes; i+=2)
    {
      final_order[i+ind] = -init_order[offset]-1;
      final_order[i+1+ind] = -init_order[offset]-1;
      offset += 2;
    }
    ind += nodes;
    for (int i = 0; i < nodes; i++)
      final_order[i+ind] = (i%2) ? 0 : 1;
    ind += nodes;
    for (int i = 0; i < nodes; i++)
      final_order[i+ind] = 1;
  }
  else if (r == 1 && f)  // CASE 5
  {
    int ind = 0;
    int offset = 0;
    for (int i = 0; i < nodes; i+=2)
    {
      final_order[i] = -init_order[offset]-1;
      final_order[i+1] = -init_order[offset]-1;
      offset += 2;
    }
    ind += nodes;
    offset = 1;
    for (int i = 0 ; i < nodes; i+=2)
    {
      final_order[i+ind] = init_order[offset];
      final_order[i+1+ind] = init_order[offset];
      offset += 2;
    }
    ind += nodes;
    for (int i = 0; i < nodes; i++)
      final_order[i+ind] = 1;
    ind += nodes;
    for (int i = 0; i < nodes; i++)
      final_order[i+ind] = (i%2) ? 1 : 0;
  }
  else if (r == 2 && f)  // CASE 6
  {
    int ind = 0;
    int offset = 0;
    for (int i = 0; i < nodes; i+=2)
    {
      final_order[i] = init_order[offset];
      final_order[i+1] = init_order[offset];
      offset += 2;
    }
    ind += nodes;
    offset = 1;
    for (int i = 0 ; i < nodes; i+=2)
    {
      final_order[i+ind] = init_order[offset];
      final_order[i+1+ind] = init_order[offset];
      offset += 2;
    }
    ind += nodes;
    for (int i = 0; i < nodes; i++)
      final_order[i+ind] = (i%2) ? 1 : 0;
    ind += nodes;
    for (int i = 0; i < nodes; i++)
      final_order[i+ind] = (i%2) ? 0 : 1;
  }
}

static unsigned const MAX_ND_ORDER = 10;

// This is all nodes, including the nodes associated with bounding edges
static inline int countTriNodes(int P)
{
  // each node on an edge has 1 dof
  // each node on a face has 2 dofs
  // each term in the following can be understood as follows
  // e*dofs*nodes
  // e      := # of entities of that dimension
  // dofs   := # of dofs associated with entities of that dimension
  // nodes  := # of nodes associated with entities of that dimension
  return 3*1*P + 1*2*P*(P-1)/2;
}
// This is all nodes, including the nodes associated with bounding edges and faces
static inline int countTetNodes(int P)
{
  // each node on an edge has 1 dof
  // each node on a face has 2 dofs
  // each node on a tet has 3 dofs
  // each term in the following can be understood as follows
  // e*dofs*nodes
  // e      := # of entities of that dimension
  // dofs   := # of dofs associated with entities of that dimension
  // nodes  := # of nodes associated with entities of that dimension
  if (P<2) return 6*1*P + 4*2*P*(P-1)/2;
  return 6*1*P + 4*2*P*(P-1)/2 + 1*3*P*(P-1)*(P-2)/6;
}

static void computeTriangleTi(
    int P, /*order*/
    mth::Matrix<double>& Q, /*Q in QR factorization of Ti*/
    mth::Matrix<double>& R) /*R in QR factorization of Ti*/
{
  const double tk[8] = {1.,0.,  -1.,1.,  0.,-1.,  0.,1.};
  const double c = 1./3.;
  int non = countTriNodes(P);

  apf::NewArray<double> eop;
  apf::NewArray<double> iop;
  getOpenPoints(P - 1, eop);
  if (P > 1)
    getOpenPoints(P - 2, iop);


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

  // Face loops to get nodes and dof2tk for faces
  for (int j = 0; j <= pm2; j++) {
    for (int i = 0; i + j <= pm2; i++)
      {
        double w = iop[i] + iop[j] + iop[pm2-i-j];
        nodes[o][0] = iop[i]/w;  nodes[o][1] = iop[j]/w;  nodes[o][2] = 0.;
        dof2tk[o++] = 0;
        nodes[o][0] = iop[i]/w;  nodes[o][1] = iop[j]/w;  nodes[o][2] = 0.;
        dof2tk[o++] = 3;
      }
  }

  // Populate T
  mth::Matrix<double> T(non,non);
  for (int m = 0; m < non; m++)
  {
    const double *tm = tk + 2*dof2tk[m];
    o = 0;

    double x = nodes[m][0]; double y = nodes[m][1];

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
  mth::decomposeQR(T, Q, R);
}

static void computeTetTi(
    int P, /*order*/
    mth::Matrix<double>& Q, /*Q in QR factorization of Ti*/
    mth::Matrix<double>& R) /*R in QR factorization of Ti*/
{
  int non = countTetNodes(P);
  const double c = 1./4.;
  const double tk[21] = { /* edge directions in a tet */
     1. ,  0.,  0.,
    -1. ,  1.,  0.,
     0. , -1.,  0.,
     0. ,  0.,  1.,
    -1. ,  0.,  1.,
     0.,  -1.,  1.,
     // this last one is the negative of the third one
     // and it is used for face and tet tangents only
     0.,   1.,  0.};


  apf::NewArray<double> eop;
  apf::NewArray<double> fop;
  apf::NewArray<double> iop;
  getOpenPoints(P - 1, eop);
  if (P > 1)
    getOpenPoints(P - 2, fop);
  if (P > 2)
    getOpenPoints(P - 3, iop);

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
  for (int i = 0; i < p; i++)  // (0,1)
  {
    nodes[o][0] = eop[i];  nodes[o][1] = 0.;  nodes[o][2] = 0.;
    dof2tk[o++] = 0;
  }
  for (int i = 0; i < p; i++)  // (1,2)
  {
    nodes[o][0] = eop[pm1-i];  nodes[o][1] = eop[i];  nodes[o][2] = 0.;
    dof2tk[o++] = 1;
  }
  for (int i = 0; i < p; i++)  // (2,0)
  {
    nodes[o][0] = 0.;  nodes[o][1] = eop[pm1-i];  nodes[o][2] = 0.;
    dof2tk[o++] = 2;
  }
  for (int i = 0; i < p; i++)  // (0,3)
  {
    nodes[o][0] = 0.;  nodes[o][1] = 0.;  nodes[o][2] = eop[i];
    dof2tk[o++] = 3;
  }
  for (int i = 0; i < p; i++)  // (1,3)
  {
    nodes[o][0] = eop[pm1-i];  nodes[o][1] = 0.;  nodes[o][2] = eop[i];
    dof2tk[o++] = 4;
  }
  for (int i = 0; i < p; i++)  // (2,3)
  {
    nodes[o][0] = 0.;  nodes[o][1] = eop[pm1-i];  nodes[o][2] = eop[i];
    dof2tk[o++] = 5;
  }

  // Face loops to get nodes and dof2tk for faces
  // (0,1,2)
  for (int j = 0; j <= pm2; j++) {
    for (int i = 0; i + j <= pm2; i++)
    {
      double w = fop[i] + fop[j] + fop[pm2-i-j];
      nodes[o][0] = fop[i]/w;  nodes[o][1] = fop[j]/w;  nodes[o][2] = 0.;
      dof2tk[o++] = 0;
      nodes[o][0] = fop[i]/w;  nodes[o][1] = fop[j]/w;  nodes[o][2] = 0.;
      dof2tk[o++] = 6;
    }
  }
  // (0,1,3)
  for (int j = 0; j <= pm2; j++) {
    for (int i = 0; i + j <= pm2; i++)
    {
      double w = fop[i] + fop[j] + fop[pm2-i-j];
      nodes[o][0] = fop[i]/w;  nodes[o][1] = 0.;  nodes[o][2] = fop[j]/w;
      dof2tk[o++] = 0;
      nodes[o][0] = fop[i]/w;  nodes[o][1] = 0.;  nodes[o][2] = fop[j]/w;
      dof2tk[o++] = 3;
    }
  }
  // (1,2,3)
  for (int j = 0; j <= pm2; j++) {
    for (int i = 0; i + j <= pm2; i++)
    {
      double w = fop[i] + fop[j] + fop[pm2-i-j];
      nodes[o][0] = fop[pm2-i-j]/w;  nodes[o][1] = fop[i]/w;  nodes[o][2] = fop[j]/w;
      dof2tk[o++] = 1;
      nodes[o][0] = fop[pm2-i-j]/w;  nodes[o][1] = fop[i]/w;  nodes[o][2] = fop[j]/w;
      dof2tk[o++] = 4;
    }
  }
  // (0,2,3)
  for (int j = 0; j <= pm2; j++) {
    for (int i = 0; i + j <= pm2; i++)
    {
      double w = fop[i] + fop[j] + fop[pm2-i-j];
      nodes[o][0] = 0.;  nodes[o][1] = fop[i]/w;  nodes[o][2] = fop[j]/w;
      dof2tk[o++] = 6;
      nodes[o][0] = 0.;  nodes[o][1] = fop[i]/w;  nodes[o][2] = fop[j]/w;
      dof2tk[o++] = 3;
    }
  }

  // Region loops to get nodes and dof2tk for regions
  for (int k = 0; k <= pm3; k++) {
    for (int j = 0; j + k <= pm3; j++) {
      for (int i = 0; i + j + k <= pm3; i++) {
        double w = iop[i] + iop[j] + iop[k] + iop[pm3-i-j-k];
        nodes[o][0] = iop[i]/w;  nodes[o][1] = iop[j]/w;  nodes[o][2] = iop[k]/w;
        dof2tk[o++] = 0;
        nodes[o][0] = iop[i]/w;  nodes[o][1] = iop[j]/w;  nodes[o][2] = iop[k]/w;
        dof2tk[o++] = 6;
        nodes[o][0] = iop[i]/w;  nodes[o][1] = iop[j]/w;  nodes[o][2] = iop[k]/w;
        dof2tk[o++] = 3;
      }
    }
  }

  // Populate T
  mth::Matrix<double> T(non,non);
  for (int m = 0; m < non; m++)
  {
    const double *tm = tk + 3*dof2tk[m];
    o = 0;

    double x = nodes[m][0]; double y = nodes[m][1]; double z = nodes[m][2];

    getChebyshevT(pm1, x, &shape_x[0]);
    getChebyshevT(pm1, y, &shape_y[0]);
    getChebyshevT(pm1, z, &shape_z[0]);
    getChebyshevT(pm1, 1. - x - y - z, &shape_l[0]);

    for (int k = 0; k <= pm1; k++) {
      for (int j = 0; j + k <= pm1; j++) {
        for (int i = 0; i + j + k <= pm1; i++)
        {
          double s = shape_x[i]*shape_y[j]*shape_z[k]*shape_l[pm1-i-j-k];
          T(o++, m) = s * tm[0];
          T(o++, m) = s * tm[1];
          T(o++, m) = s * tm[2];
        }
      }
    }
    for (int k = 0; k <= pm1; k++) {
      for (int j = 0; j + k <= pm1; j++) {
        {
          double s = shape_x[pm1-j-k]*shape_y[j]*shape_z[k];
          T(o++, m) = s*((y - c)*tm[0] - (x - c)*tm[1]);
          T(o++, m) = s*((z - c)*tm[0] - (x - c)*tm[2]);
        }
      }
    }
    for (int k = 0; k <= pm1; k++)
    {
      T(o++, m) =
        shape_y[pm1-k]*shape_z[k]*((z - c)*tm[1] - (y - c)*tm[2]);
    }
  }

  mth::decomposeQR(T, Q, R);
}

static void getTi(
    int P,
    int type,
    mth::Matrix<double>& Q,
    mth::Matrix<double>& R)
{

  bool cond = (type == apf::Mesh::TRIANGLE || type == apf::Mesh::TET);
  PCU_ALWAYS_ASSERT_VERBOSE(cond,
      "type should be either apf::Mesh::TRIANGLE or apf::Mesh::TET!");

  static apf::NewArray<double> transformQ[apf::Mesh::TYPES][MAX_ND_ORDER+1];
  static apf::NewArray<double> transformR[apf::Mesh::TYPES][MAX_ND_ORDER+1];
  int n = type == apf::Mesh::TRIANGLE ? countTriNodes(P) : countTetNodes(P);

  // get the transform matrices if the are not already computed
  if (!transformQ[type][P].allocated()) {
    mth::Matrix<double> LQ(n,n);
    mth::Matrix<double> LR(n,n);
    type == apf::Mesh::TRIANGLE ?
    	  computeTriangleTi(P, LQ, LR) : computeTetTi(P, LQ, LR);

    transformQ[type][P].allocate(n*n);
    transformR[type][P].allocate(n*n);

    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
	      transformQ[type][P][i*n+j] = LQ(i,j);
      	transformR[type][P][i*n+j] = LR(i,j);
      }
    }
  }

  // set Q and R using transformQ and transformR
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Q(i,j) = transformQ[type][P][i*n+j];
      R(i,j) = transformR[type][P][i*n+j];
    }
  }
}

template<int P>
class Nedelec: public FieldShape {
  public:
    Nedelec()
    {
      std::stringstream ss;
      ss << "Nedelec_" << P;
      name = ss.str();
      registerSelf(name.c_str());
    }
    const char* getName() const { return name.c_str(); }
    bool isVectorShape() {return true;}
    class Vertex : public apf::EntityShape
    {
    public:
      void getValues(apf::Mesh*, apf::MeshEntity*,
	  apf::Vector3 const&, apf::NewArray<double>& values) const
      {
	(void)values;
      }
      void getLocalGradients(apf::Mesh*, apf::MeshEntity*,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
        PCU_ALWAYS_ASSERT_VERBOSE(0, "error: Nedelec shapes are not \
            implemented for vertices. Aborting()!");
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
    private:
      const int dim = 1; // ref elem dim
      const double c = 0.; // center of edge
    public:
      int getOrder() {return P;}
      void getValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<double>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getValues not implemented for \
      	    Nedelec Edges. Aborting()!");
      }
      void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getLocalGradients not \
      	    implemented for Nedelec Edges. Aborting()!");
      }
      int countNodes() const {return P;}
      void alignSharedNodes(apf::Mesh*,
	  apf::MeshEntity*, apf::MeshEntity*, int order[])
      {
	(void)order;
      }
      void getVectorValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	// : to be completed
      }
      void getLocalVectorCurls(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	// : to be completed
      }
    };
    class Triangle : public apf::EntityShape
    {
    private:
      const int dim = 2; // reference element dim
      const double c = 1./3.; // center of tri
    public:
      int getOrder() {return P;}
      void getValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<double>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getValues not implemented \
      	    for Nedelec Triangle. Try getVectorValues. Aborting()!");
      }
      void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getLocalGradients not \
      	    implemented for Nedelec Triangle. Aborting()!");
      }
      int countNodes() const {return countTriNodes(P);}
      void alignSharedNodes(apf::Mesh* m,
	  apf::MeshEntity* elem, apf::MeshEntity* shared, int order[])
      {
	      int which, rotate;
	      bool flip;
        getAlignment(m,elem,shared,which,flip,rotate);
        if(!flip)
          for(int i = 0; i < P; ++i)
            order[i] = i;
        else
          for(int i = 0; i < P; ++i)
            order[i] = -(P-1-i)-1;
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


	      mth::Matrix<double> Q(dof, dof);
      	mth::Matrix<double> R(dof, dof);
      	getTi(P, apf::Mesh::TRIANGLE, Q, R);

        mth::Matrix<double> S(dof, dim);
      	for(int i = 0; i < dim; i++) // S = Ti * u
        {
          mth::Vector<double> B (dof);
          mth::Vector<double> X (dof);
          for (int j = 0; j < dof; j++) B[j] = u(j,i); // populate b in QR x = b
          mth::solveFromQR(Q, R, B, X);
          for (int j = 0; j < dof; j++)  S(j,i) = X[j]; // populate S with x
        }

        shapes.allocate(dof);
        for (int i = 0; i < dof; i++) // populate y
        {
      	  shapes[i] = apf::Vector3( S(i,0), S(i,1), 0.0 );
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
        mth::Vector<double> curlu(dof);

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
                          shape_x[i]*dshape_l[l]) * shape_y[j];
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

      	mth::Matrix<double> Q(dof, dof);
      	mth::Matrix<double> R(dof, dof);
      	getTi(P, apf::Mesh::TRIANGLE, Q, R);

        mth::Vector<double> X(dof);
        mth::solveFromQR(Q, R, curlu, X);
      }
    };
    class Tetrahedron : public apf::EntityShape
    {
    private:
      const int dim = 3;
      const double c = 1./4.;
    public:
      int getOrder() {return P;}
      void getValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<double>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getValues not implemented for \
      	    Nedelec Tetrahedron. Try getVectorValues. Aborting()!");
      }
      void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getLocalGradients not \
      	    implemented for Nedelec Tetrahedron. Aborting()!");
      }
      int countNodes() const {return countTetNodes(P);}
      void alignSharedNodes(apf::Mesh* m,
	  apf::MeshEntity* elem, apf::MeshEntity* shared, int order[])
      {
	      int which, rotate;
	      bool flip;
        getAlignment(m,elem,shared,which,flip,rotate);
        if(m->getType(shared) == apf::Mesh::EDGE)
        {
          if(!flip)
            for(int i = 0; i < P; ++i)
              order[i] = i;
          else
            for(int i = 0; i < P; ++i)
              order[i] = -(P-1-i)-1;
          return;
        }
        //must be a triangle
        int non = P*(P-1); // nodes on the face only
        int unique_non = non/2; // bc 2 dofs per node on the face
        if ( P > 1 ) // face nodes appear for k > 1
        {
          int size = (-1 + sqrt(1+8*unique_non) )/2;

          mth::Matrix<double> Nodes(size,size); // populate nodes matrix
          Nodes.zero();
          int num = 0;
          for ( int r = size-1; r >= 0; r--)
            for ( int c = size - r -1; c < size; c++)
              Nodes(r,c) = (num++);

          // CASES
          if (rotate == 1  &&  !flip) // CASE 2
          {
            for (int c = size-1; c >= 0; c--) {    // vertical swaps
              int left = size-c-1; int right = size-1;
              for(int range = left; range <= left + (right-left)/2; range++) {
                          std::swap( Nodes(range,c), Nodes(left+right-range,c) );
              }
            }
            mth::Matrix<double> temp(size, size);  // transpose
            mth::transpose(Nodes,temp);
            for (int r = 0; r < size; r++)
              for(int c = 0; c < size; c++)
                Nodes(r,c) = temp(r,c);
          }
          if (rotate == 2  &&  !flip) // CASE 3
          {
            for (int r = size-1; r >= 0; r--) {    // horizontal swaps
              int left = size-r-1; int right = size-1;
              for(int range = left; range <= left + (right-left)/2; range++) {
                          std::swap( Nodes(r,range), Nodes(r,left+right-range) );
              }
            }
            mth::Matrix<double> temp(size, size);  // transpose
            mth::transpose(Nodes,temp);
            for (int r = 0; r < size; r++)
              for(int c = 0; c < size; c++)
                Nodes(r,c) = temp(r,c);
          }
          else if (!rotate  &&  flip) // CASE 4
          {
            mth::Matrix<double> temp(size, size);  // transpose
            mth::transpose(Nodes,temp);
            for (int r = 0; r < size; r++)
              for(int c = 0; c < size; c++)
                Nodes(r,c) = temp(r,c);
          }
          else if (rotate == 1  &&  flip) // CASE 5
          {
            for (int r = size-1; r >= 0; r--) {    // horizontal swaps
              int left = size-r-1; int right = size-1;
              for(int range = left; range <= left + (right-left)/2; range++) {
                          std::swap( Nodes(r,range), Nodes(r,left+right-range) );
              }
            }
          }
          else if (rotate == 2  &&  flip) // CASE 6
          {
            for (int c = size-1; c >= 0; c--) {    // vertical swaps
              int left = size-c-1; int right = size-1;
              for(int range = left; range <= left + (right-left)/2; range++) {
                          std::swap( Nodes(range,c), Nodes(left+right-range,c) );
              }
            }
          }

          // get the init ordered list with all face nodes
          int init_order[non];
          int final_order[4*non];
          int i = 0;
          for ( int r = size-1; r >= 0; r--)
            for (int c = size-r-1 ; c < size; c++) {
              init_order[i++] = Nodes(r,c)*2;
              init_order[i++] = Nodes(r,c)*2 + 1;
            }
          alignFaceNodes(init_order, final_order, non, rotate, flip);

          for (int i = 0; i < 4*non; i++)
            order[i] = final_order[i];
        }
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

      	mth::Matrix<double> Q(dof, dof);
      	mth::Matrix<double> R(dof, dof);
      	getTi(P, apf::Mesh::TET, Q, R);

        mth::Matrix<double> S(dof, dim);
      	for(int i = 0; i < dim; i++) // S = Ti * u
        {
          mth::Vector<double> B (dof);
          mth::Vector<double> X (dof);
          for (int j = 0; j < dof; j++) B[j] = u(j,i); // populate b
          mth::solveFromQR(Q, R, B, X);
          for (int j = 0; j < dof; j++)  S(j,i) = X[j]; // populate S with x
        }

        shapes.allocate(dof);
        for (int i = 0; i < dof; i++) // populate y
        {
      	  shapes[i] = apf::Vector3( S(i,0), S(i,1), S(i,2) );
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
                      (dshape_y[j]*(y - c) + shape_y[j])*shape_x[i]*shape_z[k]);
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


      	mth::Matrix<double> Q(dof, dof);
	      mth::Matrix<double> R(dof, dof);
	      getTi(P, apf::Mesh::TET, Q, R);

        mth::Matrix<double> S(dof, dim);
      	for(int i = 0; i < dim; i++) // S = Ti * u
        {
          mth::Vector<double> B(dof);
          mth::Vector<double> X(dof);
          for (int j = 0; j < dof; j++) B[j] = u(j,i); // populate b
          mth::solveFromQR(Q, R, B, X);
          for (int j = 0; j < dof; j++)  S(j,i) = X[j]; // populate S with x
        }

        curl_shapes.allocate(dof);
        for (int i = 0; i < dof; i++) // populate y
        {
      	  curl_shapes[i] = apf::Vector3( S(i,0), S(i,1), S(i,2) );
        }
      }
    };
    EntityShape* getEntityShape(int type)
    {
      static Edge edge;
      static Triangle triangle;
      static Tetrahedron tet;
      static EntityShape* shapes[apf::Mesh::TYPES] =
      {NULL,       //vertex
       &edge,      //edge
       &triangle,  //triangle
       NULL,       //quad
       &tet,       //tetrahedron
       NULL,       //hex
       NULL,       //prism
       NULL        //pyramid
      };
      return shapes[type];
    }
    // For the following to member functions we only need to
    // consider the interior nodes, i.e.,
    // Faces: no need to count the nodes associated with bounding edges
    // Tets: no need to count the nodes associated with bounding edges/faces
    bool hasNodesIn(int dimension)
    {
      if (dimension == 1) return true;
      if (dimension == 2) return P > 1;
      if (dimension == 3) return P > 2;
      // if not returned by now dimension should be 0, so return false
      return false;
    }
    int countNodesOn(int type)
    {
      if (type == apf::Mesh::EDGE) return P;
      if (type == apf::Mesh::TRIANGLE) return 2*P*(P-1)/2;
      if (type == apf::Mesh::TET && P>2) return 3*P*(P-1)*(P-2)/6;
      // for any other case (type and P combination) return 0;
      return 0;
    }
    int getOrder() {return P;}
    void getNodeXi(int type, int node, Vector3& xi)
    {
      apf::NewArray<double> op;
      if(type == Mesh::EDGE)
      {
      	if (P > 0)
      	  getOpenPoints(P-1, op);
        xi = Vector3( -1 + 2*op[node], 0., 0. ); // map from [0,1] to [-1,1]
        return;
      }
      else if (type == Mesh::TRIANGLE)
      {
      	if (P > 1)
      	  getOpenPoints(P-2, op);
        int pm2 = P - 2; int c = 0;

        for (int j = 0; j <= pm2; j++) {
          for (int i = 0; i + j <= pm2; i++) {
            if (node/2 == c) {  // since 2 dofs per node on the face
              double w = op[i] + op[j] + op[pm2-i-j];
              xi = Vector3( op[i]/w, op[j]/w, 0. );
              return;
            }
            else
              c++;
          }
        }
      }
      else if (type == Mesh::TET)
      {
      	if (P>2)
      	  getOpenPoints(P-3, op);
        int pm3 = P - 3; int c = 0;

        for (int k = 0; k <= pm3; k++) {
          for (int j = 0; j + k <= pm3; j++) {
            for (int i = 0; i + j + k <= pm3; i++) {
              if( node/3 == c) {  // since 3 dofs per node on interior tet
                double w = op[i] + op[j] + op[k] + op[pm3-i-j-k];
                xi = Vector3( op[i]/w, op[j]/w,  op[k]/w );
                return;
              }
              else
                c++;
            }
          }
        }
      }
      else
        xi = Vector3(0,0,0);
    }
    void getNodeTangent(int type, int node, Vector3& t)
    {
      if(type == Mesh::EDGE)
      {
      	// Edges are parametrized from -1 to 1.
      	// Having the 2 here enables us to avoid multiplying dofs by 2 for
      	// when computing them for edges
        t = Vector3( 2., 0., 0.);
        return;
      }
      else if(type == Mesh::TRIANGLE)
      {
        PCU_ALWAYS_ASSERT_VERBOSE(P >= 2,
            "face nodes appear only for order bigger than or equal to 2!");
        (node % 2 == 0) ? t = Vector3(1., 0., 0.) : t = Vector3(0., 1., 0.);
        return;
      }
      else if(type == Mesh::TET)
      {
        PCU_ALWAYS_ASSERT_VERBOSE(P >= 3,
            "volume nodes appear only for order bigger than or equal to 3!");
       if (node % 3 == 0) t = Vector3(1., 0., 0.);
       else if (node % 3 == 1) t = Vector3(0., 1., 0.);
       else t = Vector3(0., 0., 1.);
       return;
      }
      else
        t = Vector3(0, 0, 0);
    }
  private:
    std::string name;
};

apf::FieldShape* getNedelec(int order)
{
  PCU_ALWAYS_ASSERT_VERBOSE(order >= 1,
      "order is expected to be bigger than or equal to 1!");
  PCU_ALWAYS_ASSERT_VERBOSE(order <= 10,
      "order is expected to be less than or equal to 10!");
  static Nedelec<1>  ND1;
  static Nedelec<2>  ND2;
  static Nedelec<3>  ND3;
  static Nedelec<4>  ND4;
  static Nedelec<5>  ND5;
  static Nedelec<6>  ND6;
  static Nedelec<7>  ND7;
  static Nedelec<8>  ND8;
  static Nedelec<9>  ND9;
  static Nedelec<10> ND10;
  static FieldShape* const nedelecShapes[10] =
  {&ND1, &ND2, &ND3, &ND4, &ND5, &ND6, &ND7, &ND8, &ND9, &ND10};
  return nedelecShapes[order-1];
}

void projectNedelecField(Field* to, Field* from)
{
  // checks on the from field
  // checks on the to field
  apf::FieldShape* tShape = getShape(to);
  std::string      tName  = tShape->getName();
  int              tOrder = tShape->getOrder();
  PCU_ALWAYS_ASSERT_VERBOSE((tName == std::string("Linear")) && (tOrder == 1),
  		"The to field needs to be 1st order Lagrange!");

  Mesh* m = getMesh(from);
  // auxiliary count fields
  Field* count = createField(m, "counter", SCALAR, getLagrange(1));
  double xis[4][3] = {{0., 0., 0.},
  	              {1., 0., 0.},
  	              {0., 1., 0.},
  	              {0., 0., 1.}};
  // zero out the fields
  zeroField(to);
  zeroField(count);

  MeshEntity* e;
  MeshIterator* it = m->begin(m->getDimension());
  while( (e = m->iterate(it)) ) {
    MeshElement* me = createMeshElement(m, e);
    Element* el = createElement(from, me);
    MeshEntity* dvs[4];
    m->getDownward(e, 0, dvs);
    for (int i=0; i<4; i++) {
      Vector3 atXi;
      getVector(el, Vector3(xis[i]), atXi);
      Vector3 currentVal;
      getVector(to, dvs[i], 0, currentVal);
      double currentCount = getScalar(count, dvs[i], 0);
      currentVal += atXi;
      currentCount += 1.;
      setVector(to, dvs[i], 0, currentVal);
      setScalar(count, dvs[i], 0, currentCount);
    }
    destroyElement(el);
    destroyMeshElement(me);
  }
  m->end(it);

  // take care of entities on part boundary
  accumulate(to);
  accumulate(count);

  it = m->begin(0);
  while( (e = m->iterate(it)) ) {
    Vector3 sum;
    getVector(to, e, 0, sum);
    setVector(to, e, 0, sum/getScalar(count, e, 0));
  }
  m->end(it);

  // take care of entities on part boundary
  synchronize(to);

  m->removeField(count);
  destroyField(count);
}

};

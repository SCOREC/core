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

// This is used for static tables only.
// The following implementation are for general orders but template classes
// are only instantiated up to order 10 for now.
static unsigned const MAX_ORDER = 10;

static inline int countEdgeNodes(int P)
{
  return (P+1);
}

static inline int countTriNodes(int P)
{
  return (P+1)*(P+2)/2;
}

static inline int countTetNodes(int P)
{
  return (P+1)*(P+2)*(P+3)/6;
}

static inline int countInternalEdgeNodes(int P)
{
  return (P-1);
}

static inline int countInternalTriNodes(int P)
{
  return (P-1)*(P-2)/2;
}

static inline int countInternalTetNodes(int P)
{
  return (P-1)*(P-2)*(P-3)/6;
}

static void computeTriangleTi(
    int P, /*order*/
    mth::Matrix<double>& Q, /*Q in QR factorization of Ti*/
    mth::Matrix<double>& R) /*R in QR factorization of Ti*/
{
  int non = countTriNodes(P);

  apf::NewArray<double> cp;
  getClosedPoints(P, cp);


  const int p = P;
  apf::NewArray<double> shape_x(p+1);
  apf::NewArray<double> shape_y(p+1);
  apf::NewArray<double> shape_l(p+1);

  apf::DynamicArray<apf::Vector3> nodes (non);

  int o = 0;

  // vertices
  nodes[o][0] = cp[0]; nodes[o][1] = cp[0]; nodes[o][2] = 0.;
  o++;
  nodes[o][0] = cp[p]; nodes[o][1] = cp[0]; nodes[o][2] = 0.;
  o++;
  nodes[o][0] = cp[0]; nodes[o][1] = cp[p]; nodes[o][2] = 0.;
  o++;

  // edges
  for (int i = 1; i < p; i++) // (0,1)
  {
    nodes[o][0] = cp[i];  nodes[o][1] = cp[0];  nodes[o][2] = 0.;
    o++;
  }

  for (int i = 1; i < p; i++) // (1,2)
  {
    nodes[o][0] = cp[p-i];  nodes[o][1] = cp[i];  nodes[o][2] = 0.;
    o++;
  }

  for (int i = 1; i < p; i++) // (2,0)
  {
    nodes[o][0] = cp[0];  nodes[o][1] = cp[p-i];  nodes[o][2] = 0.;
    o++;
  }


  // face
  for (int j = 1; j < p; j++) {
    for (int i = 1; i + j < p; i++)
      {
        double w = cp[i] + cp[j] + cp[p-i-j];
        nodes[o][0] = cp[i]/w;  nodes[o][1] = cp[j]/w;  nodes[o][2] = 0.;
        o++;
      }
  }

  // Populate T
  mth::Matrix<double> T(non,non);
  for (int m = 0; m < non; m++)
  {
    o = 0;

    double x = nodes[m][0]; double y = nodes[m][1];

    getChebyshevT(p, x, &shape_x[0]);
    getChebyshevT(p, y, &shape_y[0]);
    getChebyshevT(p, 1. - x - y, &shape_l[0]);

    for (int j = 0; j <= p; j++)
      for (int i = 0; i + j <= p; i++)
      	T(o++, m) = shape_x[i]*shape_y[j]*shape_l[p-i-j];
  }
  mth::decomposeQR(T, Q, R);
}

static void computeTetTi(
    int P, /*order*/
    mth::Matrix<double>& Q, /*Q in QR factorization of Ti*/
    mth::Matrix<double>& R) /*R in QR factorization of Ti*/
{
  int non = countTetNodes(P);

  apf::NewArray<double> cp;
  getClosedPoints(P, cp);

  const int p = P;

  apf::NewArray<double> shape_x(p+1);
  apf::NewArray<double> shape_y(p+1);
  apf::NewArray<double> shape_z(p+1);
  apf::NewArray<double> shape_l(p+1);

  apf::DynamicArray<apf::Vector3> nodes(non);
  nodes.setSize(non);


  int o = 0;
  // vertices
  nodes[o][0] = cp[0]; nodes[o][1] = cp[0]; nodes[o][2] = cp[0];
  o++;
  nodes[o][0] = cp[p]; nodes[o][1] = cp[0]; nodes[o][2] = cp[0];
  o++;
  nodes[o][0] = cp[0]; nodes[o][1] = cp[p]; nodes[o][2] = cp[0];
  o++;
  nodes[o][0] = cp[0]; nodes[o][1] = cp[0]; nodes[o][2] = cp[p];
  o++;

  // edges
  for (int i = 1; i < p; i++) // (0,1)
  {
    nodes[o][0] = cp[i];  nodes[o][1] = cp[0];  nodes[o][2] = cp[0];
    o++;
  }

  for (int i = 1; i < p; i++) // (1,2)
  {
    nodes[o][0] = cp[p-i];  nodes[o][1] = cp[i];  nodes[o][2] = cp[0];
    o++;
  }

  for (int i = 1; i < p; i++) // (2,0)
  {
    nodes[o][0] = cp[0];  nodes[o][1] = cp[p-i];  nodes[o][2] = cp[0];
    o++;
  }

  for (int i = 1; i < p; i++) // (0,3)
  {
    nodes[o][0] = cp[0];  nodes[o][1] = cp[0];  nodes[o][2] = cp[i];
    o++;
  }

  for (int i = 1; i < p; i++) // (1,3)
  {
    nodes[o][0] = cp[p-i];  nodes[o][1] = cp[0];  nodes[o][2] = cp[i];
    o++;
  }

  for (int i = 1; i < p; i++) // (2,3)
  {
    nodes[o][0] = cp[0];  nodes[o][1] = cp[p-i];  nodes[o][2] = cp[i];
    o++;
  }


  // faces
  // (0,1,2)
  for (int j = 1; j < p; j++)
    for (int i = 1; i + j < p; i++) {
      double w = cp[i] + cp[j] + cp[p-i-j];
      nodes[o][0] = cp[i]/w;  nodes[o][1] = cp[j]/w;  nodes[o][2] = cp[0];
      o++;
    }

  // (0,1,3)
  for (int j = 1; j < p; j++)
    for (int i = 1; i + j < p; i++) {
      double w = cp[i] + cp[j] + cp[p-i-j];
      nodes[o][0] = cp[i]/w;  nodes[o][1] = cp[0]/w;  nodes[o][2] = cp[j]/w;
      o++;
    }

  // (1,2,3)
  for (int j = 1; j < p; j++)
    for (int i = 1; i + j < p; i++) {
      double w = cp[i] + cp[j] + cp[p-i-j];
      nodes[o][0] = cp[p-i-j]/w;  nodes[o][1] = cp[i]/w;  nodes[o][2] = cp[j]/w;
      o++;
    }

  // (0,2,3)
  for (int j = 1; j < p; j++)
    for (int i = 1; i + j < p; i++) {
      double w = cp[i] + cp[j] + cp[p-i-j];
      nodes[o][0] = cp[0];  nodes[o][1] = cp[i]/w;  nodes[o][2] = cp[j]/w;
      o++;
    }

  // Region
  for (int k = 1; k < p; k++) {
    for (int j = 1; j + k < p; j++) {
      for (int i = 1; i + j + k < p; i++) {
        double w = cp[i] + cp[j] + cp[k] + cp[p-i-j-k];
        nodes[o][0] = cp[i]/w;  nodes[o][1] = cp[j]/w;  nodes[o][2] = cp[k]/w;
        o++;
      }
    }
  }

  // Populate T
  mth::Matrix<double> T(non,non);
  for (int m = 0; m < non; m++)
  {
    o = 0;
    double x = nodes[m][0]; double y = nodes[m][1]; double z = nodes[m][2];

    getChebyshevT(p, x, &shape_x[0]);
    getChebyshevT(p, y, &shape_y[0]);
    getChebyshevT(p, z, &shape_z[0]);
    getChebyshevT(p, 1. - x - y - z, &shape_l[0]);

    for (int k = 0; k <= p; k++)
      for (int j = 0; j + k <= p; j++)
        for (int i = 0; i + j + k <= p; i++)
          T(o++, m) = shape_x[i]*shape_y[j]*shape_z[k]*shape_l[p-i-j-k];
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

  static apf::NewArray<double> transformQ[apf::Mesh::TYPES][MAX_ORDER+1];
  static apf::NewArray<double> transformR[apf::Mesh::TYPES][MAX_ORDER+1];
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

// internal nodes only
static apf::Vector3 getH1NodeXi(int type, int P, int node)
{
  if (type == apf::Mesh::VERTEX)
    return apf::Vector3(0., 0., 0.);

  apf::NewArray<double> cp;
  getClosedPoints(P, cp);

  if (type == apf::Mesh::EDGE) {
    PCU_ALWAYS_ASSERT(node >= 0 && node < countInternalEdgeNodes(P));
    int c = 0;
    for (int i = 1; i < P; i++)
      if (node == c)
      	return apf::Vector3(2*cp[i]-1, 0., 0.);
      else
      	c++;
  }

  if (type == apf::Mesh::TRIANGLE) {
    PCU_ALWAYS_ASSERT(node >= 0 && node < countInternalTriNodes(P));
    int c = 0;
    for (int j = 1; j < P; j++)
      for (int i = 1; i + j < P; i++)
	if (node == c) {
	  double w = cp[i] + cp[j] + cp[P-i-j];
	  return apf::Vector3(cp[i]/w, cp[j]/w, 0.);
	}
	else
	  c++;
  }

  if (type == apf::Mesh::TET) {
    PCU_ALWAYS_ASSERT(node >= 0 && node < countInternalTetNodes(P));
    int c = 0;
    for (int k = 1; k < P; k++)
      for (int j = 1; j + k < P; j++)
	for (int i = 1; i + j + k < P; i++)
	  if (c == node) {
	    double w = cp[i] + cp[j] + cp[k] + cp[P-i-j-k];
	    return apf::Vector3(cp[i]/w, cp[j]/w, cp[k]/w);
	  }
	  else
	    c++;
  }

  PCU_ALWAYS_ASSERT_VERBOSE(0, "Unsupported type!");
  return apf::Vector3(0., 0., 0.);
}

template<int P>
class H1Shape: public FieldShape {
  public:
    H1Shape()
    {
      std::stringstream ss;
      ss << "H1Shape_" << P;
      name = ss.str();
      registerSelf(name.c_str());
    }
    const char* getName() const { return name.c_str(); }
    bool isVectorShape() {return false;}
    class Vertex : public apf::EntityShape
    {
    public:
      int getOrder() {return P;}
      void getValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const& /*xi*/, apf::NewArray<double>& shapes) const
      {
        shapes.allocate(1);
        shapes[0] = 1.;
      }
      void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getLocalGradients not \
      	    implemented for H1Shape for Verts. Aborting()!");
      }
      int countNodes() const {return 1;}
      void getVectorValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getVectorValues not implemented \
      	    for H1Shape. Try getValues. Aborting()!");
      }
    };
    class Edge : public apf::EntityShape
    {
    public:
      int getOrder() {return P;}
      void getValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const& xi, apf::NewArray<double>& shapes) const
      {
      	const int p = P;
      	apf::NewArray<double> shape_x(p+1);
        int dof = countNodes();

        double x = (xi[0]+1.)/2.; // go from [-1,1] to [0,1]

        poly1dBasisBarycentric(p, x, &shape_x[0]);
        shapes.allocate(dof);
        shapes[0] = shape_x[0];
        shapes[1] = shape_x[p];
        for (int i = 1; i < p; i++)
          shapes[i+1] = shape_x[i];
      }
      void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getLocalGradients not \
      	    implemented for H1Shape for Edges. Aborting()!");
      }
      int countNodes() const {return countEdgeNodes(P);}
      void getVectorValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getVectorValues not implemented \
      	    for H1Shape. Try getValues. Aborting()!");
      }
    };
    class Triangle : public apf::EntityShape
    {
    public:
      int getOrder() {return P;}
      void getValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const& xi, apf::NewArray<double>& shapes) const
      {
        const int p = P;

      	apf::NewArray<double> shape_x(p+1);
      	apf::NewArray<double> shape_y(p+1);
      	apf::NewArray<double> shape_l(p+1);

        int dof = countNodes();
        mth::Matrix<double> u(dof, 1);

        double x = xi[0]; double y = xi[1];

        getChebyshevT(p, x, &shape_x[0]);
        getChebyshevT(p, y, &shape_y[0]);
        getChebyshevT(p, 1. - x - y, &shape_l[0]);

        int n = 0;
        for (int j = 0; j <= p; j++)
          for (int i = 0; i + j <= p; i++)
            u(n++, 0) = shape_x[i]*shape_y[j]*shape_l[p-i-j];

	      mth::Matrix<double> Q(dof, dof);
      	mth::Matrix<double> R(dof, dof);
      	getTi(P, apf::Mesh::TRIANGLE, Q, R);

        mth::Matrix<double> S(dof, 1);
      	for(int i = 0; i < 1; i++) // S = Ti * u
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
      	  shapes[i] = S(i,0);
        }
      }
      void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getLocalGradients not \
      	    implemented for H1Shape for Tris. Aborting()!");
      }
      int countNodes() const {return countTriNodes(P);}
      void alignSharedNodes(apf::Mesh* m,
      	  apf::MeshEntity* elem, apf::MeshEntity* shared, int order[])
      {
      	int which, rotate;
      	bool flip;
      	getAlignment(m, elem, shared, which, flip, rotate);
      	if (!flip)
      	  for (int i = 0; i < P-1; i++)
      	    order[i] = i;
	else
	  for (int i = 0; i < P-1; i++)
	    order[i] = P-2-i;
      }
      void getVectorValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getVectorValues not implemented \
      	    for H1Shape. Try getValues. Aborting()!");
      }
    };
    class Tetrahedron : public apf::EntityShape
    {
    public:
      int getOrder() {return P;}
      void getValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const& xi, apf::NewArray<double>& shapes) const
      {
        const int p = P;

      	apf::NewArray<double> shape_x(p+1);
      	apf::NewArray<double> shape_y(p+1);
      	apf::NewArray<double> shape_z(p+1);
      	apf::NewArray<double> shape_l(p+1);

        int dof = countNodes();
        mth::Matrix<double> u(dof, 1);

        double x = xi[0]; double y = xi[1]; double z = xi[2];

        getChebyshevT(p, x, &shape_x[0]);
        getChebyshevT(p, y, &shape_y[0]);
        getChebyshevT(p, z, &shape_z[0]);
        getChebyshevT(p, 1. - x - y - z, &shape_l[0]);

        int n = 0;
        for (int k = 0; k <= p; k++)
          for (int j = 0; j + k <= p; j++)
            for (int i = 0; i + j + k <= p; i++)
            	u(n++, 0) = shape_x[i]*shape_y[j]*shape_z[k]*shape_l[p-i-j-k];


      	mth::Matrix<double> Q(dof, dof);
      	mth::Matrix<double> R(dof, dof);
      	getTi(P, apf::Mesh::TET, Q, R);

        mth::Matrix<double> S(dof, 1);
      	for(int i = 0; i < 1; i++) // S = Ti * u
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
      	  shapes[i] = S(i,0);
        }
      }
      void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getLocalGradients not \
      	    implemented for H1Shape for Tets. Aborting()!");
      }
      int countNodes() const {return countTetNodes(P);}
      void alignSharedNodes(apf::Mesh* m,
      	  apf::MeshEntity* elem, apf::MeshEntity* shared, int order[])
      {
      	int stype = m->getType(shared);
      	int which, rotate;
      	bool flip;
      	getAlignment(m, elem, shared, which, flip, rotate);
      	if (stype == apf::Mesh::EDGE) {
	  if (!flip)
	    for (int i = 0; i < P-1; i++)
	      order[i] = i;
	  else
	    for (int i = 0; i < P-1; i++)
	      order[i] = P-2-i;
	  return;
	}
	PCU_ALWAYS_ASSERT_VERBOSE(stype == apf::Mesh::TRIANGLE,
	    "shared type must be triangle!");
	int idx0, idx1;
	if (!flip) {
	  idx0 = (3-rotate) % 3;
	  idx1 = (4-rotate) % 3;
	}
	else {
	  idx0 = (2+rotate) % 3;
	  idx1 = (1+rotate) % 3;
	}
	int idx = 0;
	for (int i = 0; i <= P-3; i++)
	  for (int j = 0; j <= P-3-i; j++) {
	    int ijk[3] = {i, j, P-3-i-j};
	    order[idx] = ijk[idx0]*(P-2)-ijk[idx0]*(ijk[idx0]-1)/2+ijk[idx1];
	    idx++;
	  }
      }
      void getVectorValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
	  apf::Vector3 const& /*xi*/, apf::NewArray<apf::Vector3>& /*shapes*/) const
      {
      	PCU_ALWAYS_ASSERT_VERBOSE(0, "error: getVectorValues not implemented for \
      	    H1Shape. Try getValues. Aborting()!");
      }
    };

    EntityShape* getEntityShape(int type)
    {
      static Vertex vert;
      static Edge edge;
      static Triangle tri;
      static Tetrahedron tet;
      static apf::EntityShape* shapes[apf::Mesh::TYPES] =
      {&vert,
       &edge,
       &tri,
       NULL,
       &tet,
       NULL,
       NULL,
       NULL};
      return shapes[type];
    }
    bool hasNodesIn(int dimension)
    {
      return P > dimension;
      /* switch (dimension) { */
	/* case 0: */
	  /* return true; */
	/* case 1: */
	  /* return P>1; */
	/* case 2: */
	  /* return P>2; */
	/* case 3; */
	  /* return P>3 */
	/* default: */
	  /* return false; */
      /* } */
    }
    int countNodesOn(int type)
    {
      switch (type) {
	case apf::Mesh::VERTEX:
	  return 1;
	case apf::Mesh::EDGE:
	  if (P>1) return countInternalEdgeNodes(P); else return 0;
	case apf::Mesh::TRIANGLE:
	  if (P>2) return countInternalTriNodes(P); else return 0;
	case apf::Mesh::TET:
	  if (P>3) return countInternalTetNodes(P); else return 0;
	default:
	  return 0;
      }
    }
    int getOrder() {return P;}
    void getNodeXi(int type, int node, Vector3& xi)
    {
      xi = getH1NodeXi(type, P, node);
    }
  private:
    std::string name;
};

apf::FieldShape* getH1Shape(int order)
{
  PCU_ALWAYS_ASSERT_VERBOSE(order > 0,
      "order is expected to be bigger than or equal to 1!");
  PCU_ALWAYS_ASSERT_VERBOSE(order <= 10,
      "order is expected to be less than or equal to 10!");
  // Note: to have higher order H1 fields all you need to do is to
  // instantiate the class  up to that order in the following table
  // and change the above assert so that the code does not fail.
  static H1Shape<1>  h1_1;
  static H1Shape<2>  h1_2;
  static H1Shape<3>  h1_3;
  static H1Shape<4>  h1_4;
  static H1Shape<5>  h1_5;
  static H1Shape<6>  h1_6;
  static H1Shape<7>  h1_7;
  static H1Shape<8>  h1_8;
  static H1Shape<9>  h1_9;
  static H1Shape<10> h1_10;
  static FieldShape* const h1Shapes[11] = {NULL,
                                           &h1_1,
                                           &h1_2,
                                           &h1_3,
                                           &h1_4,
                                           &h1_5,
                                           &h1_6,
                                           &h1_7,
                                           &h1_8,
                                           &h1_9,
                                           &h1_10};
  return h1Shapes[order];
}

/* void projectL2Field(Field* to, Field* from) */
/* { */
/*   // checks on the from field */
/*   // checks on the to field */
/*   apf::FieldShape* tShape = getShape(to); */
/*   std::string      tName  = tShape->getName(); */
/*   int              tOrder = tShape->getOrder(); */
/*   PCU_ALWAYS_ASSERT_VERBOSE((tName == std::string("Linear")) && (tOrder == 1), */
/*                 "The to field needs to be 1st order Lagrange!"); */

/*   Mesh* m = getMesh(from); */
/*   // auxiliary count fields */
/*   Field* count = createField(m, "counter", SCALAR, getLagrange(1)); */
/*   double xis[4][3] = {{0., 0., 0.}, */
/*                       {1., 0., 0.}, */
/*                       {0., 1., 0.}, */
/*                       {0., 0., 1.}}; */
/*   // zero out the fields */
/*   zeroField(to); */
/*   zeroField(count); */

/*   int nc = countComponents(to); */
/*   NewArray<double> atXi(nc); */
/*   NewArray<double> currentVal(nc); */
/*   NewArray<double> sum(nc); */

/*   MeshEntity* e; */
/*   MeshIterator* it = m->begin(m->getDimension()); */
/*   while( (e = m->iterate(it)) ) { */
/*     MeshElement* me = createMeshElement(m, e); */
/*     Element* el = createElement(from, me); */
/*     MeshEntity* dvs[4]; */
/*     m->getDownward(e, 0, dvs); */
/*     for (int i=0; i<4; i++) { */
/*       getComponents(el, Vector3(xis[i]), &(atXi[0])); */
/*       getComponents(to, dvs[i], 0, &(currentVal[0])); */
/*       for (int j = 0; j < nc; j++) { */
/*         currentVal[j] += atXi[j]; */
/*       } */
/*       double currentCount = getScalar(count, dvs[i], 0); */
/*       currentCount += 1.; */
/*       setComponents(to, dvs[i], 0, &(currentVal[0])); */
/*       setScalar(count, dvs[i], 0, currentCount); */
/*     } */
/*     destroyElement(el); */
/*     destroyMeshElement(me); */
/*   } */
/*   m->end(it); */

/*   // take care of entities on part boundary */
/*   accumulate(to); */
/*   accumulate(count); */

/*   it = m->begin(0); */
/*   while( (e = m->iterate(it)) ) { */
/*     getComponents(to, e, 0, &(sum[0])); */
/*     int cnt = getScalar(count, e, 0); */
/*     for (int i = 0; i < nc; i++) { */
/*       sum[i] /= cnt; */
/*     } */
/*     setComponents(to, e, 0, &(sum[0])); */
/*   } */
/*   m->end(it); */

/*   // take care of entities on part boundary */
/*   synchronize(to); */

/*   m->removeField(count); */
/*   destroyField(count); */
/* } */

}
